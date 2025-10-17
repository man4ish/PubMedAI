import os
import json
import faiss
import numpy as np
import ollama
from drug_target_kg import add_structured_data_to_kg  # Import your KG module

# ---------------------------
# Paths
# ---------------------------
BASE_PATH = "/Volumes/Seagate2TB/PubMed"
INDEX_FOLDER = os.path.join(BASE_PATH, "Index")
INDEX_FILE = os.path.join(INDEX_FOLDER, "pubmed_index.faiss")
ID_MAP_FILE = os.path.join(INDEX_FOLDER, "pmid_map.json")
ABSTRACT_FOLDER = os.path.join(BASE_PATH, "Abstracts")

# ---------------------------
# Load FAISS index and mapping
# ---------------------------
index = faiss.read_index(INDEX_FILE)
with open(ID_MAP_FILE, "r") as f:
    pmid_list = json.load(f)

# ---------------------------
# Ollama model
# ---------------------------
MODEL_NAME = "deepseek-r1:latest"

# ---------------------------
# Functions
# ---------------------------
def get_query_embedding(query):
    result = ollama.embeddings(model=MODEL_NAME, prompt=query)
    return np.array(result["embedding"], dtype="float32").reshape(1, -1)

def get_abstract_text(pmid):
    with open(os.path.join(ABSTRACT_FOLDER, f"{pmid}.json"), "r", encoding="utf-8") as f:
        data = json.load(f)
        return data.get("abstract", "")

def extract_structured_drug_info(drug_name, abstract):
    """
    Uses DeepSeek to extract structured drug-target info from abstract.
    Returns a list of target dictionaries.
    """
    prompt = f"Extract structured drug-target-cancer information from this abstract for drug '{drug_name}':\n{abstract}\nFormat as JSON with 'target', 'cancer', and 'mechanism'."
    response = ollama.chat(model=MODEL_NAME, messages=[{"role":"user","content":prompt}])
    raw_text = response.message if hasattr(response, "message") else str(response)
    
    # Attempt to parse JSON from DeepSeek output
    try:
        targets = json.loads(raw_text)
        if isinstance(targets, dict):
            targets = [targets]  # wrap single dict as list
    except Exception:
        # fallback: return empty
        targets = []
    return targets

# ---------------------------
# User input
# ---------------------------
drug_name = input("Enter drug name: ")
query_embedding = get_query_embedding(drug_name)

# ---------------------------
# Search local index
# ---------------------------
k = 10  # top 10 relevant abstracts
faiss.normalize_L2(query_embedding)
D, I = index.search(query_embedding, k)

top_pmids = [pmid_list[i] for i in I[0]]
print(f"Found {len(top_pmids)} relevant abstracts for '{drug_name}'")

# ---------------------------
# Extract and structure drug-target info
# ---------------------------
structured_data = []

for pmid in top_pmids:
    abstract = get_abstract_text(pmid)
    targets = extract_structured_drug_info(drug_name, abstract)
    
    structured_data.append({
        "drug": drug_name,
        "targets": targets,
        "pmid": pmid,
        "abstract": abstract
    })

# ---------------------------
# Add structured data to Neo4j KG
# ---------------------------
add_structured_data_to_kg(structured_data)

print("Data successfully added to the knowledge graph!")
