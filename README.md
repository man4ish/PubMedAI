
# PubMedAI — RAG-Enabled Drug Knowledge Graph

This project builds a Retrieval-Augmented Generation (RAG) pipeline over PubMed abstracts to explore drug–target–disease relationships and visualize them as a knowledge graph.

## Project Structure

```
PubMedAI/
├── config.py                      # Global configuration (paths, constants, etc.)
├── README.md                      # Project documentation
│
├── data_pipeline/                 # Data ingestion & indexing
│   ├── get_pubmed_abstract.py     # Fetches abstracts from PubMed
│   └── index_abstracts.py         # Builds FAISS or text index for retrieval
│
├── retriever/                     # Embedding & retrieval logic
│   ├── embedding_utils.py         # Embedding generator (e.g., using sentence-transformers)
│   └── faiss_retriever.py         # Handles vector search with FAISS
│
├── generator/                     # LLM or RAG-based text generation
│   └── rag_generator.py           # Generates answers using retrieved documents
│
├── kg/                            # Knowledge graph utilities
│   └── drug_target_kg.py          # Builds or visualizes a drug-target knowledge graph
│
└── scripts/                       # Entry-point scripts for CLI use
    ├── query_pubmed.py            # Main interactive query script
    └── visualize_drug_kg.py       # Visualizes drug-target graph
```

## Installation

1. Clone the repository
   ```bash
   git clone https://github.com/<your-username>/PubMedAI.git
   cd PubMedAI
   ```

2. Create and activate environment
   ```bash
   conda create -n pubmedai python=3.10 -y
   conda activate pubmedai
   ```

3. Install dependencies
   ```bash
   pip install -r requirements.txt
   ```

## Workflow Overview

### Step 1 — Download PubMed Abstracts
Fetch abstracts related to biomedical keywords (e.g., drugs, diseases, or proteins).
```bash
python data_pipeline/get_pubmed_abstract.py
```
Output: Saves `.txt` files under `data/Abstracts/`.

### Step 2 — Create FAISS Index
Convert abstracts into embeddings and index them using FAISS for fast semantic search.
```bash
python data_pipeline/index_abstracts.py
```
Output: Stores vector index (`pubmed_index.faiss`) and metadata (`pmid_map.json`) in `data/Index/`.

### Step 3 — Query the Knowledge Base
Query by drug name, mechanism, or disease — uses semantic retrieval (FAISS) to find related abstracts.
```bash
python scripts/query_pubmed.py
```
Example:
```text
Enter drug name or question: tyrosine kinase inhibitor
```
Output:
- Lists relevant abstracts
- Updates `data/knowledge_graph.json`

### Step 4 — Visualize the Knowledge Graph
Build and visualize a graph of drug–target–disease relationships.
```bash
python scripts/visualize_drug_kg.py
```
Output: Interactive graph (nodes = drugs, targets, diseases).

## How It Works (RAG Architecture)

| Component                | Description                                                                  |
|-------------------------|------------------------------------------------------------------------------|
| Retriever               | FAISS vector store indexes PubMed embeddings for efficient semantic search   |
| Augmenter               | Retrieved abstracts enrich context for querying and relationship extraction  |
| Generator (optional)    | Can be extended to use an LLM (e.g., GPT) for natural-language summarization |
| Knowledge Graph         | Represents entities (drugs, genes, diseases) and their inferred connections  |

## Example Flow
```
Input: "tyrosine kinase inhibitor"
↓
FAISS → retrieves 10 top abstracts
↓
NER → extracts entities (Drug, Target, Cancer)
↓
Graph → adds new nodes/edges
↓
Visualization → interactive network
```

## Future Enhancements
- Integrate OpenAI or Llama for summarization (RAG generation stage)
- Add biomedical NER (e.g., SciSpacy or BERN2)
- Build web dashboard for interactive KG exploration
- Add temporal trends and citation analytics
```