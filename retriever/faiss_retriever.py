import faiss
import json
import numpy as np
import os
from config import INDEX_FILE, ID_MAP_FILE

class FAISSRetriever:
    def __init__(self):
        self.index = faiss.read_index(INDEX_FILE)
        with open(ID_MAP_FILE, "r") as f:
            self.pmid_list = json.load(f)

    def search(self, query_embedding, top_k=10):
        faiss.normalize_L2(query_embedding)
        D, I = self.index.search(query_embedding, top_k)
        return [self.pmid_list[i] for i in I[0]]
