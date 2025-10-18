import os
import json
import faiss
import numpy as np
import ollama

# ---------------------------
# Paths
# ---------------------------
BASE_PATH = "/Volumes/Seagate2TB/PubMed"
ABSTRACT_FOLDER = os.path.join(BASE_PATH, "Abstracts")
INDEX_FOLDER = os.path.join(BASE_PATH, "Index")
INDEX_FILE = os.path.join(INDEX_FOLDER, "pubmed_index.faiss")
ID_MAP_FILE = os.path.join(INDEX_FOLDER, "pmid_map.json")

os.makedirs(INDEX_FOLDER, exist_ok=True)

# ---------------------------
# Load abstracts
# ---------------------------
abstracts = []
pmid_list = []

for file in os.listdir(ABSTRACT_FOLDER):
    if file.endswith(".json"):
        with open(os.path.join(ABSTRACT_FOLDER, file), "r", encoding="utf-8") as f:
            data = json.load(f)
            text = data.get("abstract", "")
            if text.strip():
                abstracts.append(text)
                pmid_list.append(data["pmid"])

print(f"Loaded {len(abstracts)} abstracts for indexing.")

# ---------------------------
# Generate embeddings using Ollama
# ---------------------------
model_name = "deepseek-r1:latest"
embeddings = []
batch_size = 20  # adjust if you hit memory issues

for i in range(0, len(abstracts), batch_size):
    batch = abstracts[i:i+batch_size]
    for abstract in batch:
        result = ollama.embeddings(model=model_name, prompt=abstract)
        embeddings.append(result["embedding"])
    print(f"Processed {i + len(batch)} / {len(abstracts)} abstracts")

embeddings = np.array(embeddings, dtype="float32")

# ---------------------------
# Normalize for cosine similarity
# ---------------------------
faiss.normalize_L2(embeddings)

# ---------------------------
# Build FAISS index
# ---------------------------
dimension = embeddings.shape[1]
index = faiss.IndexFlatIP(dimension)
index.add(embeddings)

print(f"FAISS index built with {index.ntotal} vectors.")

# ---------------------------
# Save index and ID map
# ---------------------------
faiss.write_index(index, INDEX_FILE)
with open(ID_MAP_FILE, "w", encoding="utf-8") as f:
    json.dump(pmid_list, f, indent=2)

print(f"Index saved to: {INDEX_FILE}")
print(f"PMID map saved to: {ID_MAP_FILE}")
