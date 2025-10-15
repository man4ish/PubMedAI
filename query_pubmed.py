from Bio import Entrez
import ollama
import numpy as np
import faiss
import json
import os

# ---------------------------
# Paths
# ---------------------------
BASE_PATH = "/Volumes/Seagate2TB/PubMed"
ABSTRACT_FOLDER = os.path.join(BASE_PATH, "Abstracts")
INDEX_FOLDER = os.path.join(BASE_PATH, "Index")
INDEX_FILE = os.path.join(INDEX_FOLDER, "pubmed_index.faiss")
ID_MAP_FILE = os.path.join(INDEX_FOLDER, "pmid_map.json")

# ---------------------------
# Load FAISS index and ID map
# ---------------------------
index = faiss.read_index(INDEX_FILE)
with open(ID_MAP_FILE, "r", encoding="utf-8") as f:
    pmid_list = json.load(f)

# ---------------------------
# Ollama client
# ---------------------------
model_name = "deepseek-r1:latest"

# ---------------------------
# Search local FAISS index
# ---------------------------
def search_local(query, top_k=5):
    query_embedding = ollama.embeddings(model=model_name, prompt=query)["embedding"]
    query_embedding = np.array([query_embedding], dtype="float32")
    faiss.normalize_L2(query_embedding)
    scores, indices = index.search(query_embedding, top_k)
    
    results = []
    for idx in indices[0]:
        pmid = pmid_list[idx]
        file_path = os.path.join(ABSTRACT_FOLDER, f"{pmid}.json")
        if os.path.exists(file_path):
            with open(file_path, "r", encoding="utf-8") as f:
                data = json.load(f)
                results.append(data.get("abstract", ""))
    return results

# ---------------------------
# Search PubMed online
# ---------------------------
def search_pubmed_online(query, max_results=5, email="your_email@example.com"):
    Entrez.email = email
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    pmids = record["IdList"]
    
    abstracts = []
    for pmid in pmids:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        fetch_record = Entrez.read(handle)
        handle.close()
        article = fetch_record['PubmedArticle'][0]['MedlineCitation']['Article']
        abstract_text = article.get('Abstract', {}).get('AbstractText', [""])[0]
        abstracts.append(abstract_text)
    return abstracts

# ---------------------------
# RAG Answer
# ---------------------------
def rag_answer(query, top_k_local=5, top_k_online=5):
    local_docs = search_local(query, top_k_local)
    online_docs = search_pubmed_online(query, top_k_online)
    combined_docs = local_docs + online_docs
    context = "\n\n".join(combined_docs)
    
    prompt = f"""You are a biomedical research assistant.
Use the following abstracts to answer the question.

Question: {query}

Context:
{context}

Answer:"""
    
    response = ollama.generate(model=model_name, prompt=prompt)
    return response["response"]

# ---------------------------
# Example usage
# ---------------------------
query = "AI applications in cancer genomics"
answer = rag_answer(query)
print(answer)
