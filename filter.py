import os
import requests
import logging
import pandas as pd
from Bio import Entrez
import nltk
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
import string
import time  # For rate limiting
import arxiv
from nltk.stem import PorterStemmer
from Bio import Medline
from io import StringIO
import threading
from scholarly import scholarly
import xml.etree.ElementTree as ET
from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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
        logging.info(f"Searching PubMed for: {query}")
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results,
            usehistory="y"  # Use history to fetch records later
        )
        results = Entrez.read(handle)
        handle.close()
        logging.info(f"Number of articles found: {results['Count']}")
        return results
    except Exception as e:
        logging.error(f"An error occurred during PubMed search: {e}")
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
            articles.extend(list(Medline.parse(StringIO(data))))
            logging.info(f"Fetching records {start + 1} to {end}")
        return articles
    except Exception as e:
        logging.error(f"An error occurred while fetching PubMed details: {e}")
        return []

def fetch_all_arxiv_articles(query, total_results=208483, batch_size=100):
    """
    Fetch all articles from arXiv based on the query with enhanced logging and progress tracking.
    
    Args:
        query (str): The search query.
        total_results (int): Total number of results to fetch.
        batch_size (int): Number of results per request.
        
    Returns:
        list: A list of dictionaries containing article details.
    """
    arxiv_articles = []
    start = 0
    retries = 3  # Number of retries for failed requests

    with tqdm(total=total_results, desc="Fetching arXiv Articles") as pbar:
        while start < total_results:
            attempt = 0
            success = False
            while attempt < retries and not success:
                try:
                    logging.info(f"Requesting page: Start={start}, Max Results={batch_size}")
                    search_query = f"{query} start={start} max_results={batch_size}"
                    search = arxiv.Search(
                        query=search_query,
                        max_results=batch_size,
                        sort_by=arxiv.SortCriterion.Relevance
                    )

                    batch = []
                    for result in search.results():
                        arxiv_article = {
                            'Title': result.title,
                            'Authors': ', '.join([author.name for author in result.authors]),
                            'Abstract': result.summary,
                            'PublicationDate': result.published.strftime('%Y-%m-%d'),
                            'DOI': result.doi if result.doi else '',
                            'Source': 'arXiv'
                        }
                        batch.append(arxiv_article)

                    if not batch:
                        logging.info("No more articles found. Ending fetch.")
                        success = True
                        break

                    arxiv_articles.extend(batch)
                    pbar.update(len(batch))
                    logging.info(f"Fetched {len(batch)} articles. Total so far: {len(arxiv_articles)}")

                    # Increment start for next batch
                    start += batch_size

                    # Sleep to respect arXiv's rate limits (e.g., 3 seconds between requests)
                    sleep_duration = 3
                    logging.info(f"Sleeping for {sleep_duration} seconds to respect rate limits.")
                    time.sleep(sleep_duration)

                    success = True  # Exit retry loop

                except Exception as e:
                    attempt += 1
                    logging.error(f"Attempt {attempt} failed: {e}")
                    logging.info(f"Retrying after sleep. Attempt {attempt} of {retries}.")
                    time.sleep(5)  # Sleep before retrying

            if not success:
                logging.error(f"Failed to fetch batch starting at {start} after {retries} attempts. Skipping to next batch.")
                start += batch_size
                pbar.update(batch_size)  # Skip updating for failed batch

    return arxiv_articles

def validate_full_text(paper):
    """
    Validate if the paper meets the full-text criteria.
    
    Args:
        paper (dict): A dictionary containing paper details.
        
    Returns:
        bool: True if the paper meets the criteria, False otherwise.
    """
    # Placeholder implementation - modify based on your criteria
    required_fields = ['Title', 'Authors', 'Abstract', 'PublicationDate']
    for field in required_fields:
        if field not in paper or not paper[field]:
            return False
    return True

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
    arxiv_articles = fetch_all_arxiv_articles(search_query, total_results=208483, batch_size=100)
    if arxiv_articles:
        print(f"Number of arXiv articles retrieved: {len(arxiv_articles)}")
    else:
        arxiv_articles = []
        print("No arXiv articles found.")

    # Combine and Filter Articles
    articles = pubmed_articles + arxiv_articles
    print(f"Total number of articles retrieved: {len(articles)}")

    if articles:
        filtered_papers = filter_papers(
            articles,
            inclusion_keywords=inclusion_criteria,
            exclusion_keywords=exclusion_criteria,
            verbose=False
        )
        if filtered_papers:
            save_results_to_csv(filtered_papers)
            print(f"{len(filtered_papers)} papers match the criteria.")
        else:
            print("No papers matched the criteria.")
    else:
        print("No articles to process.")