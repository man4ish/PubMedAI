import ollama
import json
from config import MODEL_NAME, ABSTRACT_FOLDER
import os

def get_abstract_text(pmid: str) -> str:
    file_path = os.path.join(ABSTRACT_FOLDER, f"{pmid}.json")
    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return data.get("abstract", "")

def rag_answer(user_query, pmid_list):
    context = "\n\n".join([get_abstract_text(pmid) for pmid in pmid_list])
    prompt = f"Answer the following question using the context below:\nQuestion: {user_query}\n\nContext:\n{context}"
    response = ollama.chat(model=MODEL_NAME, messages=[{"role":"user","content":prompt}])
    return response.message if hasattr(response, "message") else str(response)

def extract_structured_info(drug_name, abstract):
    prompt = f"""
    Extract structured drug-target-cancer information for '{drug_name}':
    {abstract}
    Format as JSON: [{{'target': str, 'cancer': str, 'mechanism': str}}]
    """
    response = ollama.chat(model=MODEL_NAME, messages=[{"role": "user", "content": prompt}])
    text = getattr(response, "message", str(response))
    try:
        parsed = json.loads(text)
        return parsed if isinstance(parsed, list) else [parsed]
    except Exception:
        return []
