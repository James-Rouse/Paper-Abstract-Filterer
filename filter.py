from Bio import Entrez
import pandas as pd
import nltk
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
import string
import os
import time  # For rate limiting
import arxiv
from nltk.stem import PorterStemmer
from Bio import Medline
from Bio import Entrez
from io import StringIO
import threading
import re

# Define the custom nltk_data path
nltk_data_path = r'C:\Users\wooot\nltk_data'
os.environ['NLTK_DATA'] = nltk_data_path
os.makedirs(nltk_data_path, exist_ok=True)
nltk.data.path.append(nltk_data_path)

# Ensure 'punkt' tokenizer is available
try:
    nltk.data.find('tokenizers/punkt')
except LookupError:
    nltk.download('punkt', download_dir=nltk_data_path)

# Ensure 'stopwords' is available
try:
    nltk.data.find('corpora/stopwords')
except LookupError:
    nltk.download('stopwords', download_dir=nltk_data_path)

# Set your email to use Entrez (required by NCBI)
Entrez.email = "james.thomas.rouse@gmail.com"

def search_pubmed(query, max_results=20000):
    """
    Search PubMed with the given query and return search results.
    """
    try:
        print(f"Searching PubMed for: {query}")
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results,
            usehistory="y"  # Use history to fetch records later
        )
        results = Entrez.read(handle)
        handle.close()
        print(f"Number of articles found: {results['Count']}")
        return results
    except Exception as e:
        print(f"An error occurred during PubMed search: {e}")
        return None

def fetch_pubmed_details(pubmed_results):
    """
    Fetch article details from PubMed using WebEnv and QueryKey.
    """
    try:
        articles = []
        batch_size = 100
        webenv = pubmed_results["WebEnv"]
        query_key = pubmed_results["QueryKey"]
        id_list = pubmed_results["IdList"]

        for start in range(0, len(id_list), batch_size):
            end = min(start + batch_size, len(id_list))
            print(f"Fetching records {start + 1} to {end}")
            fetch_handle = Entrez.efetch(
                db="pubmed",
                rettype="medline",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=webenv,
                query_key=query_key,
            )
            data = fetch_handle.read()
            fetch_handle.close()

            # Correctly parse Medline data using StringIO
            records = list(Medline.parse(StringIO(data)))
            
            for record in records:
                article = {
                    'Title': record.get('TI', ''),
                    'Abstract': record.get('AB', ''),
                    'Authors': ', '.join(record.get('AU', [])),
                    'PublicationDate': record.get('DP', ''),
                    'DOI': record.get('LID', '').replace(' [doi]', ''),
                    'Source': 'PubMed'
                }
                articles.append(article)
            time.sleep(0.5)  # Respect rate limits
        return articles
    except Exception as e:
        print(f"An error occurred while fetching PubMed details: {e}")
        return []

def validate_full_text(paper):
    """
    Check specific details in the paper such as brain scans, lesions, and confounding factors.
    """
    abstract = paper.get('Abstract', '').lower()
    title = paper.get('Title', '').lower()
    combined_text = title + " " + abstract

    brain_scan_keywords = ["mri", "ct", "brain scan"]
    lesion_keywords = ["lesion", "clear lesion"]
    confounding_keywords = ["toxic", "drug-induced", "metabolic", "physiologic", "pharmacological"]

    # Verify presence of brain scan
    has_brain_scan = any(keyword in combined_text for keyword in brain_scan_keywords)
    # Verify presence of lesion clarity
    has_lesion = any(keyword in combined_text for keyword in lesion_keywords)
    # Exclude if confounding factors are present
    has_confounding = any(keyword in combined_text for keyword in confounding_keywords)

    return has_brain_scan and has_lesion and not has_confounding

def search_arxiv(query, max_results=500):
    """
    Search arXiv with the given query and return search results as a list of dictionaries.
    """
    try:
        print(f"Searching arXiv for: {query}")
        client = arxiv.Client(page_size=100, delay_seconds=2)  # Adjust delay_seconds and limit max_results
        search = arxiv.Search(
            query=query,
            max_results=max_results,
            sort_by=arxiv.SortCriterion.Relevance
        )
        arxiv_articles = []
        for result in client.results(search):
            arxiv_article = {
                'Title': result.title,
                'Abstract': result.summary,
                'Authors': ', '.join([author.name for author in result.authors]),
                'PublicationDate': result.published.strftime('%Y-%m-%d'),
                'DOI': result.doi or '',
                'Source': 'arXiv'
            }
            arxiv_articles.append(arxiv_article)
        print(f"Number of arXiv articles found: {len(arxiv_articles)}")
        return arxiv_articles
    except Exception as e:
        print(f"An error occurred during arXiv search: {e}")
        return []

def search_arxiv_with_timeout(query, max_results=500, timeout=30):
    """
    Search arXiv with a timeout to prevent hanging indefinitely.
    """
    result = []

    def search():
        nonlocal result
        result = search_arxiv(query, max_results)

    thread = threading.Thread(target=search)
    thread.start()
    thread.join(timeout)  # Join the thread with a timeout
    if thread.is_alive():
        print("arXiv search timed out.")
        return []
    return result

def filter_papers(papers, inclusion_keywords, exclusion_keywords, verbose=False):
    """
    Filter papers based on inclusion and exclusion criteria.
    """
    filtered = []
    for paper in papers:
        try:
            if validate_full_text(paper):
                filtered.append(paper)
            elif verbose:
                print("Excluded: Paper does not meet full-text criteria")
        except Exception as e:
            print(f"Error processing paper: {e}")
    return filtered

def save_results_to_csv(filtered_papers, file_name="filtered_papers.csv"):
    """
    Save filtered papers to a CSV file.
    """
    try:
        df = pd.DataFrame(filtered_papers)
        df.to_csv(file_name, index=False)
        print(f"Results saved to {file_name}")
    except Exception as e:
        print(f"Error saving results to CSV: {e}")

if __name__ == "__main__":
    search_query = (
        "(case reports[Publication Type] OR case series[Publication Type]) AND "
        "(delusion) AND "
        "(brain scan OR MRI OR CT)"
    )
    inclusion_criteria = ["case report", "case series", "delusion", "brain scan", "mri", "ct"]
    exclusion_criteria = ["toxic", "drug-induced", "metabolic", "physiologic", "pharmacological"]

    # PubMed Search
    pubmed_results = search_pubmed(search_query)
    if pubmed_results:
        pubmed_articles = fetch_pubmed_details(pubmed_results)
        print(f"Number of PubMed articles retrieved: {len(pubmed_articles)}")
    else:
        pubmed_articles = []
        print("No PubMed articles found.")

    # arXiv Search
    arxiv_articles = search_arxiv_with_timeout("delusion brain imaging MRI CT", max_results=500, timeout=30)
    if arxiv_articles:
        print(f"Number of arXiv articles retrieved: {len(arxiv_articles)}")
    else:
        arxiv_articles = []
        print("No arXiv articles found.")

    # Combine and Filter Articles
    articles = pubmed_articles + arxiv_articles
    print(f"Total number of articles retrieved: {len(articles)}")

    if articles:
        filtered_papers = filter_papers(articles, inclusion_keywords=inclusion_criteria, exclusion_keywords=exclusion_criteria, verbose=False)
        if filtered_papers:
            save_results_to_csv(filtered_papers)
            print(f"{len(filtered_papers)} papers match the criteria.")
        else:
            print("No papers matched the criteria.")
    else:
        print("No articles to process.")
