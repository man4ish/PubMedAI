import os
import json
import requests
from bs4 import BeautifulSoup
from Bio import Entrez

# ---------------------------
# User Config
# ---------------------------
Entrez.email = "mandecent.gupta@egmail.com"  # Replace with your email
search_term = "cancer genomics"         # Your PubMed search term
retmax = 5000                              # Number of articles to fetch

# External HDD paths
BASE_PATH = "/Volumes/Seagate2TB/PubMed"
ABSTRACT_FOLDER = os.path.join(BASE_PATH, "Abstracts")
METADATA_FOLDER = os.path.join(BASE_PATH, "Metadata")
PDF_FOLDER = os.path.join(BASE_PATH, "PDFs")

# Create folders if they don't exist
os.makedirs(ABSTRACT_FOLDER, exist_ok=True)
os.makedirs(METADATA_FOLDER, exist_ok=True)
os.makedirs(PDF_FOLDER, exist_ok=True)

# ---------------------------
# Step 1: Search PubMed
# ---------------------------
handle = Entrez.esearch(db="pubmed", term=search_term, retmax=retmax)
record = Entrez.read(handle)
id_list = record["IdList"]
print(f"Found {len(id_list)} articles for '{search_term}'")

# ---------------------------
# Step 2: Download abstracts & metadata
# ---------------------------
for pmid in id_list:
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        article = records['PubmedArticle'][0]['MedlineCitation']['Article']

        # Extract info
        title = article.get('ArticleTitle', '')
        abstract = ''
        if 'Abstract' in article:
            abstract = ' '.join([t for t in article['Abstract']['AbstractText']])

        authors = []
        if 'AuthorList' in article:
            for a in article['AuthorList']:
                name = f"{a.get('LastName','')} {a.get('ForeName','')}".strip()
                if name:
                    authors.append(name)

        # Save abstract
        abstract_data = {"pmid": pmid, "title": title, "abstract": abstract}
        with open(os.path.join(ABSTRACT_FOLDER, f"{pmid}.json"), "w", encoding="utf-8") as f:
            json.dump(abstract_data, f, ensure_ascii=False, indent=2)

        # Save metadata
        metadata = {"pmid": pmid, "title": title, "authors": authors}
        with open(os.path.join(METADATA_FOLDER, f"{pmid}_meta.json"), "w", encoding="utf-8") as f:
            json.dump(metadata, f, ensure_ascii=False, indent=2)

        print(f"Saved PMID {pmid} - {title}")

        # ---------------------------
        # Step 3: Try to download open-access PDF from PMC
        # ---------------------------
        # Build PMC search URL
        pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/?term={pmid}"
        response = requests.get(pmc_url)
        soup = BeautifulSoup(response.text, 'html.parser')
        pdf_link = soup.find('a', string='PDF')
        if pdf_link:
            pdf_href = pdf_link.get('href')
            if pdf_href.startswith('/'):
                pdf_href = "https://www.ncbi.nlm.nih.gov" + pdf_href

            pdf_response = requests.get(pdf_href)
            pdf_path = os.path.join(PDF_FOLDER, f"{pmid}.pdf")
            with open(pdf_path, 'wb') as f:
                f.write(pdf_response.content)
            print(f"Downloaded PDF for PMID {pmid}")

    except Exception as e:
        print(f"Error with PMID {pmid}: {e}")
