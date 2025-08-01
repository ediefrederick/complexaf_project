import requests
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from itertools import chain
from collections import Counter
import pandas as pd
import os
from requests.exceptions import RequestException

def all_chains_over_50(pdb_id):
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        resp = requests.get(url, timeout=10)
        if not resp.ok:
            return False

        data = resp.json()
        entity_ids = data.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
        for entity_id in entity_ids:
            ent_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
            ent_resp = requests.get(ent_url, timeout=10)
            if not ent_resp.ok:
                return False
            ent_data = ent_resp.json()
            length = ent_data.get("entity_poly", {}).get("rcsb_sample_sequence_length", 0)
            if length < 50:
                return False
        return True
    except Exception:
        return False

def is_disallowed_species(species_list):
    disallowed_keywords = ["virus", "viral", "phage"]
    for species in species_list:
        if not isinstance(species, str):
            continue
        species_lower = species.lower()
        if species_lower in ["unknown", "synthetic construct"]:
            return True
        if any(keyword in species_lower for keyword in disallowed_keywords):
            return True
    return False


def get_species_per_entity(pdb_id):
    species_info = {}
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        n_entities = len(resp.json().get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", []))

        for entity_id in range(1, n_entities + 1):
            ent_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
            ent_resp = requests.get(ent_url, timeout=10)
            ent_resp.raise_for_status()
            ent_data = ent_resp.json()
            # filename = f"{pdb_id}_{entity_id}.json"
            # # Open the file in write mode and dump the dictionary to it
            # with open(filename, 'w') as f:
            #     json.dump(ent_data, f, indent=4)
            sources = ent_data.get("rcsb_entity_source_organism", [])
            species = [src.get("ncbi_scientific_name", "Unknown") for src in sources]
            species_info[entity_id] = species if species else ["Unknown"]

    except Exception:
        return {}
    return species_info

def validate_pdb(pdb_id, query_type):
    # Returns pdb_id if valid, else None
    if query_type != "query_na":
        if not all_chains_over_50(pdb_id):
            return None

    species_info = get_species_per_entity(pdb_id)
    if any(is_disallowed_species(species) for species in species_info.values()):
        return None

    return pdb_id


def fetch_pdb_ids(query, rows_per_page=100, max_pages=100):
    all_pdb_ids = set()
    start = 0

    while True:
        query['request_options']['paginate']['start'] = start
        query['request_options']['paginate']['rows'] = rows_per_page

        response = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query)
        if response.status_code != 200:
            raise RuntimeError(f"RCSB query failed: {response.status_code} {response.text}")

        results = response.json().get("result_set", [])
        if not results:
            break

        for rec in results:
            pdb_id = rec['identifier'].split('_')[0]
            all_pdb_ids.add(pdb_id)

        start += rows_per_page
        if start // rows_per_page >= max_pages:
            print("Reached max_pages limit.")
            break
        time.sleep(0.1)

    return sorted(all_pdb_ids)

def get_mutation_info(pdb_id):
    mutations = []
    try:
        # First, get the number of polymer entities
        url_main = f"https://data.rcsb.org/rest/v1/core/structure/{pdb_id}"
        main_response = requests.get(url_main).json()
        entities = main_response.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])

        for entity_id in entities:
            url_entity = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
            entity_response = requests.get(url_entity).json()

            # Check mutation field
            mutation = entity_response.get("rcsb_entity_source_organism", {}).get("mutation")
            if mutation:
                mutations.append(f"Entity {entity_id}: {mutation}")

            # Alternatively, check details field for 'mutant'
            description = entity_response.get("entity", {}).get("pdbx_description", "")
            if "mutant" in description.lower():
                mutations.append(f"Entity {entity_id}: {description}")

    except Exception as e:
        return f"Error fetching {pdb_id}: {e}"

    return mutations if mutations else ["None found"]


def get_experimental_method(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        method = data.get("exptl", [{}])[0].get("method", "Unknown")
        return {
            "pdb_id": pdb_id,
            "experimental_method": method
        }
    except Exception as e:
        print(f"Error fetching experimental method for {pdb_id}: {e}")
        return None

def get_pisa_interfaces(pdb_id, timeout=180):
    url = f"https://www.ebi.ac.uk/pdbe/api/pisa/assembly/{pdb_id.lower()}/1"
    
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code != 200:
            print(f"Request failed for {pdb_id} with status code {response.status_code}")
            return []
        
        data = response.json()
        if pdb_id.lower() not in data or "assembly" not in data[pdb_id.lower()]:
            print(f"No assembly data for {pdb_id} assembly 1")
            return []

        entry = data[pdb_id.lower()]
        assembly_data = entry["assembly"]
        formula = assembly_data.get("formula", "unknown")

        # Skip if formula is null or contains lowercase letters
        if formula is None:
            return []
        if any(char.islower() for char in formula):
            return []

        # Clean formula
        if isinstance(formula, str):
            formula = formula.replace("(", "").replace(")", "")

        composition = assembly_data.get("composition", "unknown")
        interface_count = assembly_data.get("interface_count", "unknown")

        return [{
            "assembly_id": "1",
            "formula": formula,
            "composition": composition,
            "interface_count": interface_count
        }]

    except RequestException as e:
        print(f"Request exception for {pdb_id}: {e}")
        return []
    

def get_uniprot_mapping(pdb_id, timeout=180):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"

    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code != 200:
            print(f"Failed to get UniProt mapping for {pdb_id}")
            return {}

        data = response.json()
        chain_to_uniprot = {}

        uniprot_data = data.get(pdb_id.lower(), {}).get("UniProt", {})
        for uniprot_id, info in uniprot_data.items():
            mappings = info.get("mappings", [])
            for mapping in mappings:
                chain_id = mapping.get("chain_id")
                if chain_id:
                    chain_to_uniprot[chain_id] = uniprot_id

        return chain_to_uniprot

    except RequestException as e:
        print(f"Request exception for {pdb_id}: {e}")
        return {}


# def download_fasta_from_rcsb(pdb_id, output_dir="fasta_files", max_retries=5):
#     url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
#     retries = 0
#     os.makedirs(output_dir, exist_ok=True)
#     while retries < max_retries:
#         try:
#             response = requests.get(url, timeout=10)
#             response.raise_for_status()
#             fasta_text = response.text

#             fasta_path = os.path.join(output_dir, f"{pdb_id}.fasta")
#             with open(fasta_path, "w") as f:
#                 f.write(fasta_text)
#             print(f"Saved FASTA for {pdb_id}")
#             return True

#         except requests.exceptions.HTTPError as e:
#             if response.status_code == 429:
#                 wait_time = 5 * (retries + 1)
#                 print(f"Rate limited for {pdb_id}, sleeping {wait_time}s")
#                 time.sleep(wait_time)
#                 retries += 1
#             elif response.status_code == 404:
#                 print(f"No FASTA found for {pdb_id} (404)")
#                 return False
#             else:
#                 print(f"HTTP error for {pdb_id}: {e}")
#                 return False
#         except Exception as e:
#             print(f"Error downloading FASTA for {pdb_id}: {e}")
#             return False

#     print(f"Max retries exceeded for {pdb_id}")
#     return False


# def download_fasta_batch(pdb_ids):
#     for pdb_id in pdb_ids:
#         download_fasta_from_rcsb(pdb_id)

def main():
    query_protein = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_accession_info.deposit_date", "operator": "greater_or_equal", "value": "2023-01-30"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.diffrn_resolution_high.value", "operator": "less_or_equal", "value": 3.5}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.selected_polymer_entity_types", "operator": "exact_match", "value": "Protein (only)"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.nonpolymer_entity_count", "operator": "equals", "value": 0}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity.pdbx_description", "operator": "exact_match", "negation": True, "value": "ANTIBODY"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity.pdbx_description", "operator": "exact_match", "negation": True, "value": "NANOBODY"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name", "operator": "exact_match", "negation": True, "value": "Viruses"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.polymer_entity_count_protein", "operator": "greater_or_equal", "value": 2}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "entity_poly.rcsb_sample_sequence_length", "operator": "greater", "value": 50}},
            ]
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": 100},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "scoring_strategy": "combined"
        }
    }

    query_na = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_accession_info.deposit_date", "operator": "greater_or_equal", "value": "2023-01-30"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.diffrn_resolution_high.value", "operator": "less_or_equal", "value": 3.5}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.selected_polymer_entity_types", "operator": "exact_match", "value": "Protein/NA"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.nonpolymer_entity_count", "operator": "equals", "value": 0}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity.pdbx_description", "operator": "exact_match", "negation": True, "value": "ANTIBODY"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.polymer_entity_count_protein", "operator": "greater_or_equal", "value": 2}},
            ]
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": 100},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "scoring_strategy": "sequence"
        }
    }
    # query_type = "query_protein"
    # pdb_ids = fetch_pdb_ids(query_protein)

    # print(f"Initial results count: {len(pdb_ids)}")

    # valid_pdb_ids = []
    valid_pdb_ids = ['8CIF', '8CMR', '8CQZ', '8CRA', '8FZZ', '8G0K', '8G0P', '8G0Q', '8G19', '8G1B', '8G1C', '8G25', '8GBX', '8GC1', '8GDO', '8GJP', '8GKJ', '8GKL', '8GLC', '8IC9', '8IF0', '8IGD', '8IGU', '8ILX', '8IM0', '8IPC', '8IVW', '8J09', '8J1I', '8J4A', '8J64', '8J6K', '8J7E', '8J8P', '8J8Q', '8J9B', '8JBJ', '8JDI', '8JEL', '8JEN', '8JEO', '8JEP', '8JEQ', '8JI9', '8JJ6', '8JMQ', '8JTK', '8JTL', '8JZT', '8OI2', '8ONW', '8OWS', '8OZ3', '8OZB', '8P6H', '8P8O', '8P9D', '8PA0', '8PFD', '8PNL', '8PU3', '8PWU', '8PWV', '8PWZ', '8Q5D', '8QKR', '8QP6', '8QUS', '8R4N', '8R61', '8R7A', '8R7D', '8REK', '8RGI', '8RO5', '8RO6', '8RO7', '8RO8', '8ROC', '8RPP', '8RTX', '8RW9', '8RWF', '8RY2', '8RZ2', '8S0N', '8S2M', '8S5T', '8S83', '8SG3', '8SI1', '8SJ4', '8SKK', '8SOW', '8SOZ', '8SVG', '8SVH', '8SXP', '8T2D', '8T9O', '8T9P', '8TE7', '8TFN', '8TFR', '8TGE', '8TJF', '8TJT', '8TKT', '8TKU', '8TS7', '8TS8', '8TVJ', '8TZN', '8TZW', '8U3S', '8U70', '8U8C', '8UHT', '8UKH', '8UKI', '8ULH', '8UZP', '8V4I', '8V52', '8VBK', '8VBL', '8VBM', '8VBN', '8VBP', '8VBQ', '8VBR', '8VBS', '8VEG', '8VFU', '8VGE', '8VGF', '8VGG', '8VK1', '8VK2', '8VRS', '8VS8', '8VTE', '8VU1', '8VUA', '8VVB', '8W6B', '8W70', '8W90', '8W9J', '8WHM', '8WT4', '8XAA', '8XEH', '8XIE', '8XJ0', '8XLA', '8XTA', '8XUP', '8Y2K', '8Y31', '8Y77', '8Y9S', '8Y9T', '8Y9U', '8YBL', '8YBM', '8YBN', '8YBO', '8YBP', '8YD8', '8YDD', '8YJ3', '8YKT', '8YM5', '8YM6', '8YNX', '8YOR', '8YX1', '8YXV', '8Z22', '8Z4C', '8Z50', '8Z6G', '8Z8M', '8Z8V', '8ZD6', '8ZH7', '8ZQG', '9AVO', '9AWE', '9B7D', '9B9F', '9BDI', '9BIV', '9BIY', '9BJG', '9BJH', '9CMC', '9CQI', '9CQJ', '9D65', '9D6J', '9D6M', '9DBO', '9DDT', '9DH2', '9DL0', '9DMF', '9DQ3', '9DQ5', '9DSC', '9DX6', '9E69', '9EG6', '9EGN', '9EKD', '9EKE', '9ERU', '9ERW', '9ETL', '9EYH', '9F18', '9F5T', '9F98', '9FVC', '9FZD', '9G1Y', '9G9N', '9G9O', '9G9P', '9GNV', '9GNX', '9GP2', '9H4R', '9H4S', '9H6A', '9H6B', '9H6C', '9H6D', '9HI7', '9HIZ', '9I8A', '9II9', '9IMU', '9INC', '9IP6', '9IV3', '9JBQ', '9JWA', '9K3E', '9KO8', '9L5U', '9LHJ', '9MJ6', '9MJC', '9MJD', '9MJI', '9MK4', '9MOJ', '9MPB', '9R2Z', '9RMM', '9UX0']
    # max_workers = 20  # tweak this number depending on your bandwidth and server limits
    # validate = partial(validate_pdb, query_type=query_type)

    # with ThreadPoolExecutor(max_workers=max_workers) as executor:
    #     futures = {executor.submit(validate, pdb_id): pdb_id for pdb_id in pdb_ids}
    #     for future in as_completed(futures):
    #         pdb_id = future.result()
    #         if pdb_id:
    #             valid_pdb_ids.append(pdb_id)


    # for pdb_id in pdb_ids:
    #     species_dict = get_species_per_entity(pdb_id)
    #     flat_species_list = [sp for species_list in species_dict.values() for sp in species_list]
    #     disallowed = is_disallowed_species(flat_species_list)
    #     sequence_len = all_chains_over_50(pdb_id)
    #     if disallowed == False and sequence_len == True:
    #         valid_pdb_ids.append(pdb_id)
  
    # print(f"Validated PDB count: {valid_pdb_ids}")


    # records = []

    # for pdb_id in valid_pdb_ids:
    #     mutations_per_entity = {}
    #     mutation_list = get_mutation_info(pdb_id)
    #     if isinstance(mutation_list, list) and mutation_list != ["None found"]:
    #         for m in mutation_list:
    #             if m.startswith("Entity "):
    #                 try:
    #                     entity_id, mut = m.replace("Entity ", "").split(":", 1)
    #                     mutations_per_entity[int(entity_id.strip())] = mut.strip()
    #                 except ValueError:
    #                     continue

    #     experiment_info = get_experimental_method(pdb_id)
    #     method = experiment_info["experimental_method"] if experiment_info else "N/A"

    #     species_dict = get_species_per_entity(pdb_id)
    #     for entity_id in species_dict:
    #         entity_name = f"{pdb_id}_{entity_id}"
    #         mutation = mutations_per_entity.get(entity_id, "N/A")
    #         records.append({
    #             "PDB ID": pdb_id,
    #             "Entity ID": entity_name,
    #             "Mutations": mutation,
    #             "Experimental Method": method
    #         })

    # df = pd.DataFrame(records)
    # print(df.head())
    # df.to_csv("pdb_entity_info.csv", index=False)

    # Load existing dataframe
    df = pd.read_csv("pdb_entity_info.csv")

    # Add new column for FASTA path
    fasta_dir = "query_fasta"
    df["FASTA Path"] = df["PDB ID"].apply(
        lambda pdb_id: os.path.join(fasta_dir, f"{pdb_id}.fasta")
        if os.path.exists(os.path.join(fasta_dir, f"{pdb_id}.fasta")) else "N/A"
    )

    # Save the updated CSV
    df.to_csv("pdb_entity_info.csv", index=False)
    print("Added FASTA Path column to CSV.")


    # # Initialize new columns
    # df["Protein Species"] = "N/A"
    # df["UniProt Mapping"] = "N/A"

    # # Create a lookup from (pdb_id, entity_id) → species / uniprot
    # for pdb_id in df["PDB ID"].unique():
    #     try:
    #         # Species info is a dict: {1: [species1, species2], 2: [...], ...}
    #         species_dict = get_species_per_entity(pdb_id)
    #         uniprot_dict = get_uniprot_mapping(pdb_id)

    #         for entity_id, species_list in species_dict.items():
    #             full_entity_id = f"{pdb_id}_{entity_id}"
    #             species_str = "; ".join(species_list) if species_list else "N/A"

    #             # Filter dataframe rows for this entity
    #             mask = df["Entity ID"] == full_entity_id
    #             df.loc[mask, "Protein Species"] = species_str

    #             # Find UniProt mapping for any chain in this entity
    #             matching_chains = [chain for chain in uniprot_dict if chain.startswith(str(entity_id))]
    #             uniprot_ids = set([uniprot_dict[chain] for chain in matching_chains])
    #             if uniprot_ids:
    #                 df.loc[mask, "UniProt Mapping"] = "; ".join(sorted(uniprot_ids))

    #     except Exception as e:
    #         print(f"Skipping {pdb_id} due to error: {e}")

    # Save updated dataframe
    df.to_csv("pdb_entity_info.csv", index=False)


    # for pdb_id in valid_pdb_ids:
    #     pisa_assemblies = get_pisa_interfaces(pdb_id)
    #     uniprot = get_uniprot_mapping(pdb_id)
    #     for assembly in pisa_assemblies:
    #         formula = assembly.get("formula", "unknown")
            # print(f"{pdb_id}: {formula}")

    # download_fasta_batch(valid_pdb_ids)

if __name__ == "__main__":
    main()

