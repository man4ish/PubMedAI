import os

BASE_PATH = "/Volumes/Seagate2TB/PubMed"
ABSTRACT_FOLDER = os.path.join(BASE_PATH, "Abstracts")
INDEX_FOLDER = os.path.join(BASE_PATH, "Index")
INDEX_FILE = os.path.join(INDEX_FOLDER, "pubmed_index.faiss")
ID_MAP_FILE = os.path.join(INDEX_FOLDER, "pmid_map.json")

MODEL_NAME = "deepseek-r1:latest"
TOP_K = 10
