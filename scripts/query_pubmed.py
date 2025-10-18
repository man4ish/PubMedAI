from retriever.embedding_utils import get_query_embedding
from retriever.faiss_retriever import FAISSRetriever
from generator.rag_generator import rag_answer, get_abstract_text, extract_structured_info
from kg.drug_target_kg import add_structured_data_to_kg
from config import TOP_K

def main():
    user_query = input("Enter drug name or question: ")

    retriever = FAISSRetriever()
    query_embedding = get_query_embedding(user_query)
    top_pmids = retriever.search(query_embedding, TOP_K)
    print(f"\nFound {len(top_pmids)} relevant abstracts for '{user_query}'")

    structured_data = []
    for pmid in top_pmids:
        abstract = get_abstract_text(pmid)
        info = extract_structured_info(user_query, abstract)
        structured_data.append({"drug": user_query, "targets": info, "pmid": pmid})

    add_structured_data_to_kg(structured_data)
    print("Structured drug-target info added to the knowledge graph!\n")

    answer = rag_answer(user_query, top_pmids)
    print("=== RAG Answer ===")
    print(answer)

if __name__ == "__main__":
    main()
