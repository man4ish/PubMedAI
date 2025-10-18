import numpy as np
import ollama
from config import MODEL_NAME

def get_query_embedding(query: str) -> np.ndarray:
    result = ollama.embeddings(model=MODEL_NAME, prompt=query)
    return np.array(result["embedding"], dtype="float32").reshape(1, -1)
