from typing import Any
import httpx
from mcp.server.fastmcp import FastMCP
import asyncio
from urllib.parse import quote
 
mcp = FastMCP('cheminformatics_mcp')



valid_resources = [
    "activity", "assay", "atc_class", "binding_site", "biotherapeutic", "cell_line",
    "chembl_id_lookup", "compound_record", "compound_structural_alert", "document",
    "document_similarity", "document_term", "drug", "drug_indication", "drug_warning",
    "go_slim", "mechanism", "metabolism", "molecule", "molecule_form", "organism",
    "protein_classification", "similarity", "source", "status", "substructure", "target",
    "target_component", "target_relation", "tissue", "xref_source"
]

# Tool 1
@mcp.tool()
async def fetch_chembl_id(chembl_id: str, timeout: float = 10.0):
    """
    Fetches detailed information for a ChEMBL ID by first looking up the ID's resource type,
    then retrieving the full data from the appropriate endpoint.
    
    Parameters:
    -----------
    chembl_id : str
        A valid ChEMBL identifier
    timeout : float, optional
        Maximum time in seconds to wait for the API responses (default: 10.0)
        
    Returns:
    --------
    Optional[Dict[str, Any]]
        The complete data for the ChEMBL ID as a dictionary if successful, None otherwise
    """

    base_url = "https://www.ebi.ac.uk"
    lookup_url = f"{base_url}/chembl/api/data/chembl_id_lookup/{chembl_id}.json"
    async with httpx.AsyncClient(timeout=timeout) as client:
        try:
            lookup_response = await client.get(lookup_url)
            lookup_response.raise_for_status()
            lookup_data = lookup_response.json()

            resource_url = lookup_data.get("resource_url") + ".json"
            if not resource_url:
                print(f"No resource_url found for {chembl_id}")
                return None

            full_url = base_url + resource_url
            resource_response = await client.get(full_url)
            resource_response.raise_for_status()
            return resource_response.json()

        except httpx.HTTPStatusError as e:
            print(f"HTTP error {e.response.status_code} for {chembl_id}: {e.response.text}")
        except Exception as e:
            print(f"Unexpected error for {chembl_id}: {e}")
    return None

# Tool 2
@mcp.tool()
async def fetch_chembl_given_resources(chembl_id: str, resources: str, timeout: float = 10.0):
    """
    Fetches data directly from a specific ChEMBL API resource endpoint for a given ChEMBL ID.
    Use this function when you already know which resource type the ChEMBL ID belongs to.
    
    Parameters:
    -----------
    chembl_id : str
        A valid ChEMBL identifier
    resources : str
        The specific ChEMBL resource type to query. Must be one of the valid resources listed
        in the valid_resources list (e.g., 'molecule', 'target', 'assay')
    timeout : float, optional
        Maximum time in seconds to wait for the API response (default: 10.0)
        
    Returns:
    --------
    Optional[Dict[str, Any]]
        The data from the specified resource for the ChEMBL ID as a dictionary if successful, 
        None otherwise
    """

    if resources not in valid_resources:
        print(f"Invalid resource: {resources}. Valid resources are: {valid_resources}")
        return None

    url = f"https://www.ebi.ac.uk/chembl/api/data/{resources}/{chembl_id}.json"
    async with httpx.AsyncClient(timeout=timeout) as client:
        try:
            response = await client.get(url)
            response.raise_for_status()
            return response.json()

        except httpx.HTTPStatusError as e:
            print(f"HTTP error {e.response.status_code} for {chembl_id}: {e.response.text}")
        except Exception as e:
            print(f"Unexpected error for {chembl_id}: {e}")
    return None


# Tool 3
@mcp.tool()
async def substructure_search(molecule: str, limit: int = 10, timeout: float = 10.0, offset: int = 0):
    """
    Performs a substructure search in the ChEMBL database given a molecule.
    This function finds molecules in the ChEMBL database that contain the specified
    substructure pattern defined by the SMILES string.
    
    Parameters:
    -----------
    molecule : str
        The SMILES, CHEMBL_ID or InChIKey
    limit : int, optional
        Maximum number of results to return (default: 10)
    timeout : float, optional
        Maximum time in seconds to wait for the API response (default: 10.0)
    offset : int, optional
        Number of results to skip for pagination (default: 0)
        
    Returns:
    --------
    Optional[List[Dict[str, str]]]
        A list of dictionaries containing the ChEMBL ID and canonical SMILES 
        for each matching molecule if successful, None otherwise.
        Each dictionary has keys 'CHEMBL_ID' and 'SMILES'.
    """

    encoded_molecule = quote(molecule)
    url = f"https://www.ebi.ac.uk/chembl/api/data/substructure/{encoded_molecule}.json"
    async with httpx.AsyncClient(timeout=timeout) as client:
        try:
            response = await client.get(url, params={"limit": limit, "offset": offset})
            response.raise_for_status()

            if "application/json" in response.headers.get("Content-Type", ""):
                data = response.json()
                results = []
                for mol in data.get("molecules", []):
                    chembl_id = mol.get("molecule_chembl_id")
                    mol_structures = mol.get("molecule_structures", {})
                    canonical_smiles = mol_structures.get("canonical_smiles")
                    if chembl_id and canonical_smiles:
                        results.append({
                            "CHEMBL_ID": chembl_id,
                            "SMILES": canonical_smiles
                        })
                return results
            else:
                print("Unexpected content type:", response.text)
        except httpx.HTTPError as e:
            print(f"HTTP error: {e}")
        except Exception as e:
            print(f"Unexpected error: {e}")
    return None


# Tool 4
@mcp.tool()
async def similarity_search(molecule: str, similarity: int=80, limit: int = 10, timeout: float = 30.0, offset: int = 0):
    """
    Performs a similarity search in the ChEMBL database using a SMILES string.
    This function finds molecules in the ChEMBL database that are structurally
    similar to the specified molecule, using the provided similarity threshold.
    
    Parameters:
    -----------
    molecule : str
        The SMILES, CHEMBL ID or InChIKey
    similarity : int, optional
        Similarity threshold as a percentage (0-100)
    limit : int, optional
        Maximum number of results to return (default: 10)
    timeout : float, optional
        Maximum time in seconds to wait for the API response (default: 30.0)
    offset : int, optional
        Number of results to skip for pagination (default: 0)
        
    Returns:
    --------
    Optional[List[Dict[str, str]]]
        A list of dictionaries containing the ChEMBL ID and canonical SMILES 
        for each similar molecule if successful, None otherwise.
        Each dictionary has keys 'CHEMBL_ID' and 'SMILES'.
    """

    encoded_molecule = quote(molecule)
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{encoded_molecule}/{int(similarity)}.json"
    async with httpx.AsyncClient(timeout=timeout) as client:
        try:
            response = await client.get(url, params={"limit": limit, "offset": offset})
            response.raise_for_status()

            if "application/json" in response.headers.get("Content-Type", ""):
                data = response.json()
                results = []
                for mol in data.get("molecules", []):
                    chembl_id = mol.get("molecule_chembl_id")
                    mol_structures = mol.get("molecule_structures", {})
                    canonical_smiles = mol_structures.get("canonical_smiles")
                    if chembl_id and canonical_smiles:
                        results.append({
                            "CHEMBL_ID": chembl_id,
                            "SMILES": canonical_smiles
                        })
                return results
            else:
                print("Unexpected content type:", response.text)
        except httpx.HTTPError as e:
            print(f"HTTP error: {e}")
        except Exception as e:
            print(f"Unexpected error: {e}")
    return None

if __name__ == "__main__":
    mcp.run(transport='stdio')