# PubMedAI

**PubMedAI** is a hybrid Retrieval-Augmented Generation (RAG) pipeline for biomedical research. It enables you to:

- Download and store PubMed abstracts locally.
- Build a semantic search index using **DeepSeek embeddings**.
- Query your local index offline for fast retrieval.
- Optionally query PubMed online for new publications.
- Generate concise answers using the DeepSeek LLM based on both local and online sources.

---

## Table of Contents

1. [Project Structure](#project-structure)  
2. [Requirements](#requirements)  
3. [Quickstart](#quickstart)  
4. [Step 1: Download PubMed Abstracts](#step-1-download-pubmed-abstracts)  
5. [Step 2: Build FAISS Index](#step-2-build-faiss-index)  
6. [Step 3: Query Abstracts (RAG Style)](#step-3-query-abstracts-rag-style)  
7. [Notes](#notes)  
8. [Future Enhancements](#future-enhancements)  

---

## Project Structure

```
PubMedAI/
├── get_pubmed_abstract.py # Download PubMed abstracts
├── index_abstracts.py # Generate embeddings & build FAISS index
├── query_pubmed.py # Query local + online data using RAG
├── README.md # Project documentation
└── External HDD/
└── PubMed/
├── Abstracts/ # JSON abstracts downloaded from PubMed
└── Index/ # FAISS index + PMID mapping
```


---

## Requirements

- Python 3.10+  
- Install required packages:

```bash
pip install biopython faiss-cpu sentence-transformers ollama numpy
```

## Ollama installed with DeepSeek model:

```
ollama pull deepseek-r1:latest
```

External HDD mounted (example path: /Volumes/Seagate2TB)

## Quickstart

### Download PubMed abstracts.

- Build the FAISS index using DeepSeek embeddings.

### Query abstracts using the hybrid RAG pipeline.

- Step 1: Download PubMed Abstracts

Edit get_pubmed_abstract.py with your query term and run:

```
python get_pubmed_abstract.py
```

This will save abstracts in JSON format under:

/Volumes/Seagate2TB/PubMed/Abstracts/

- Step 2: Build FAISS Index

After downloading abstracts:

```
python index_abstracts.py
```

This script will:

Generate embeddings for all abstracts using DeepSeek.

## Build a FAISS index for fast semantic search.

Save the index and a mapping of vector → PMID to:

```
/Volumes/Seagate2TB/PubMed/Index/
```

- Step 3: Query Abstracts (RAG Style)

Run the query script:

```
python query_pubmed.py
```

### Enter your biomedical query at the prompt.

The script retrieves the most relevant abstracts from your local FAISS index.

Optionally, it can also fetch new abstracts from PubMed online.

Generates a concise answer using the DeepSeek LLM.

Example Query:

Enter your PubMed query: AI applications in cancer genomics