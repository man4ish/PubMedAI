import os
import json
import faiss
import numpy as np
import ollama

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
# Function to get embeddings for query
# ---------------------------
def get_query_embedding(query):
    result = ollama.embeddings(model=MODEL_NAME, prompt=query)
    return np.array(result["embedding"], dtype="float32").reshape(1, -1)

# ---------------------------
# Function to fetch abstract text
# ---------------------------
def get_abstract_text(pmid):
    with open(os.path.join(ABSTRACT_FOLDER, f"{pmid}.json"), "r", encoding="utf-8") as f:
        data = json.load(f)
        return data.get("abstract", "")

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
# Extract drug-target info using DeepSeek
# ---------------------------
results = []

for pmid in top_pmids:
    abstract = get_abstract_text(pmid)
    prompt = f"Extract drug-target information from this abstract for drug '{drug_name}':\n{abstract}"
    
    response = ollama.chat(model=MODEL_NAME, messages=[{"role":"user","content":prompt}])
    # DeepSeek returns 'message' field
    summary = response.message if hasattr(response, "message") else str(response)
    
    results.append({
        "pmid": pmid,
        "abstract": abstract,
        "drug": drug_name,
        "drug_target_summary": summary
    })

# ---------------------------
# Display results
# ---------------------------
for r in results:
    print("\n---")
    print(f"PMID: {r['pmid']}")
    print(f"Drug: {r['drug']}")
    print(f"Drug-Target Summary: {r['drug_target_summary']}")
