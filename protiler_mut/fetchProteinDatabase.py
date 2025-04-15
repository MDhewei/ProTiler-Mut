import os
import json
import requests
from typing import Dict

def download_and_load_jsons() -> Dict[str, dict]:
    """
    Download and load all of our annotation JSONs from Figshare.

    Returns:
        A dict with keys:
            - 'exons_dic'
            - 'domain_dic'
            - 'alpha_dic'
            - 'po_dic'
            - 'clinvar_dic'
        each mapping to the loaded JSON content.
    """
    # 1) Define your Figshare download URLs here:
    url_mapping = {
        "exons_dic": "https://figshare.com/ndownloader/files/53596592",
        "domain_dic": "https://figshare.com/ndownloader/files/53596595",
        "alpha_dic": "https://figshare.com/ndownloader/files/53596607",
        "po_dic": "https://figshare.com/ndownloader/files/53596601",
        "clinvar_dic": "https://figshare.com/ndownloader/files/53596589"
    }

    # 2) (Optional) Where to cache the files locally
    cache_dir = os.path.expanduser("~/.cache/protein_annotations")
    os.makedirs(cache_dir, exist_ok=True)

    results: Dict[str, dict] = {}

    for key, url in url_mapping.items():
        local_path = os.path.join(cache_dir, f"{key}.json")

        # Download if not already cached
        if not os.path.exists(local_path):
            resp = requests.get(url, stream=True)
            resp.raise_for_status()
            with open(local_path, "wb") as fh:
                for chunk in resp.iter_content(1 << 20):
                    fh.write(chunk)
            print(f"[downloaded] {key} → {local_path}")
        else:
            print(f"[cached]    {key} → {local_path}")

        # Load JSON into memory
        with open(local_path, "r") as fh:
            results[key] = json.load(fh)

    return results


# --- Usage ---
if __name__ == "__main__":
    dicts = download_and_load_jsons()
    exons_dic   = dicts["exons_dic"]
    domain_dic  = dicts["domain_dic"]
    alpha_dic   = dicts["alpha_dic"]
    po_dic      = dicts["po_dic"]
    clinvar_dic = dicts["clinvar_dic"]

    # quick sanity check
    print("Exons entries:",   len(exons_dic))
    print("Domains entries:", len(domain_dic))
    print("AlphaFold entries:", len(alpha_dic))
    print("Phospho entries:",  len(po_dic))
    print("ClinVar entries:",  len(clinvar_dic))
