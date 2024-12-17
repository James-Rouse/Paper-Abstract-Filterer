# filepath: /c:/Users/wooot/OneDrive/Desktop/Paper Abstract Filterer/filter.py
from Bio import Entrez
import pandas as pd
import nltk
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
import string
import os
import time  # For rate limiting
import logging

# Configure logging
logging.basicConfig(
    level=logging.WARNING,  # Set to WARNING to suppress DEBUG and INFO messages
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Define the custom nltk_data path
nltk_data_path = r'C:\Users\wooot\nltk_data'
os.environ['NLTK_DATA'] = nltk_data_path
os.makedirs(nltk_data_path, exist_ok=True)
nltk.data.path.append(nltk_data_path)

# Ensure 'punkt' tokenizer is available
try:
    nltk.data.find('tokenizers/punkt')
    print("Punkt tokenizer found!")
except LookupError:
    print("Punkt tokenizer not found, downloading...")
    nltk.download('punkt', download_dir=nltk_data_path)

# Ensure 'stopwords' is available
try:
    nltk.data.find('corpora/stopwords')
    print("Stopwords found!")
except LookupError:
    print("Stopwords not found, downloading...")
    nltk.download('stopwords', download_dir=nltk_data_path)

# Set your email to use Entrez (required by NCBI)
Entrez.email = "james.thomas.rouse@gmail.com"  # Replace with your actual email address

def search_pubmed(query, max_results=10000):
    """
    Search PubMed with the given query and return search results.
    
    Args:
        query (str): The search query.
        max_results (int): Maximum number of results to retrieve (limit to 10,000).
        
    Returns:
        dict or None: Search results if successful, else None.
    """
    try:
        print(f"Searching PubMed for: {query}")
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results,       # Limit to 10,000 results
            usehistory="y"            # Use history to fetch records later
        )
        results = Entrez.read(handle)
        handle.close()
        count = int(results.get("Count", 0))
        print(f"Number of articles found: {count}")
        
        if count > max_results:
            print(f"Note: Only the first {max_results} articles will be fetched due to API limitations.")
        
        return results
    except Exception as e:
        print(f"An error occurred during PubMed search: {e}")
        return None

def fetch_pubmed_details(search_results, verbose=False):
    """
    Fetch detailed PubMed information using WebEnv and QueryKey.
    Extract and flatten the required fields from each PubmedArticle.
    
    Args:
        search_results (dict): Results from Entrez.esearch.
        verbose (bool): If True, print sample papers for inspection.
    
    Returns:
        list: A list of dictionaries with flattened paper details.
    """
    records = []
    try:
        count = int(search_results.get("Count", 0))
        webenv = search_results.get("WebEnv", "")
        query_key = search_results.get("QueryKey", "")
        batch_size = 100

        if count == 0:
            print("No articles found.")
            return records

        # Cap the total number of records to 10,000
        max_allowed = 10000
        if count > max_allowed:
            print(f"Warning: Found {count} articles, but only {max_allowed} will be fetched.")
            count = max_allowed

        for start in range(0, count, batch_size):
            end = min(start + batch_size, count)
            print(f"\rFetching records {start+1} to {end}", end='', flush=True)
            handle = Entrez.efetch(
                db="pubmed",
                rettype="xml",
                retmode="xml",
                retstart=start,
                retmax=batch_size,
                webenv=webenv,
                query_key=query_key
            )
            result = Entrez.read(handle)
            handle.close()
            time.sleep(0.4)  # Be respectful of NCBI rate limits

            # Check for 'PubmedArticle' in the result
            if 'PubmedArticle' in result:
                for article in result['PubmedArticle']:
                    # Initialize a flat dictionary
                    flat_paper = {}

                    # Extract ArticleTitle
                    try:
                        flat_paper['ArticleTitle'] = article['MedlineCitation']['Article']['ArticleTitle']
                    except KeyError:
                        flat_paper['ArticleTitle'] = ''

                    # Extract Authors
                    try:
                        authors = article['MedlineCitation']['Article']['AuthorList']
                        author_names = []
                        for author in authors:
                            if 'LastName' in author and 'Initials' in author:
                                name = f"{author['LastName']} {author['Initials']}"
                                author_names.append(name)
                            elif 'CollectiveName' in author:
                                author_names.append(author['CollectiveName'])
                        flat_paper['AuthorList'] = '; '.join(author_names)
                    except KeyError:
                        flat_paper['AuthorList'] = ''

                    # Extract Abstract
                    try:
                        abstract_parts = article['MedlineCitation']['Article']['Abstract']['AbstractText']
                        if isinstance(abstract_parts, list):
                            abstract = ' '.join(abstract_parts)
                        else:
                            abstract = abstract_parts
                        flat_paper['Abstract'] = abstract
                    except KeyError:
                        flat_paper['Abstract'] = ''

                    # Extract Publication Date from Journal -> JournalIssue -> PubDate
                    try:
                        journal = article['MedlineCitation']['Article']['Journal']
                        journal_issue = journal.get('JournalIssue', {})
                        pub_date_info = journal_issue.get('PubDate', {})
                        year = pub_date_info.get('Year', '')
                        month = pub_date_info.get('Month', '')
                        day = pub_date_info.get('Day', '')

                        # Handle cases where only MedlineDate is available
                        if not year:
                            medline_date = pub_date_info.get('MedlineDate', '')
                            if medline_date:
                                year = medline_date.split(' ')[0]  # Extract the first part as year

                        # Format PubDate
                        pub_date = f"{year}-{month}-{day}".strip('-')  # Remove leading/trailing dashes
                        flat_paper['PubDate'] = pub_date
                    except KeyError:
                        flat_paper['PubDate'] = ''

                    records.append(flat_paper)
            else:
                print("\nNo 'PubmedArticle' found in the result.")
                print("Result keys:", result.keys())

        # Print a newline character to move the cursor to the next line after fetching is complete
        print()

        # **Print the structure of a few papers for inspection**
        if verbose:
            print("\nSample papers structure:")
            for idx, paper in enumerate(records[:5], start=1):
                print(f"\nSample Paper {idx}:")
                for key, value in paper.items():
                    print(f"{key}: {value}")

        return records
    except Exception as e:
        print(f"An error occurred while fetching PubMed details: {e}")
        return records

def validate_full_text(paper, verbose=False):
    """
    Validate if the paper meets the full-text criteria.
    
    Args:
        paper (dict): A dictionary containing paper details.
        verbose (bool): If True, print exclusion reasons.
        
    Returns:
        bool: True if the paper meets the criteria, False otherwise.
    """
    # Required fields after flattening
    required_fields = ['ArticleTitle', 'AuthorList', 'PubDate']
    missing_fields = [field for field in required_fields if not paper.get(field, '').strip()]

    if missing_fields:
        if verbose:
            print(f"Paper excluded due to missing field(s): {', '.join(missing_fields)}")
        return False
    return True

def filter_papers(papers, inclusion_keywords, exclusion_keywords, verbose=False):
    """
    Filter papers based on inclusion and exclusion criteria.
    
    Args:
        papers (list): A list of paper dictionaries.
        inclusion_keywords (list): Keywords for inclusion.
        exclusion_keywords (list): Keywords for exclusion.
        verbose (bool): If True, print summary of exclusions.
        
    Returns:
        list: A list of filtered paper dictionaries.
    """
    filtered = []
    exclusion_reasons = {}

    for paper in papers:
        try:
            # Validate required fields without printing exclusion reasons here
            valid = validate_full_text(paper, verbose=False)
            if not valid:
                reason = "Missing required fields"
                exclusion_reasons[reason] = exclusion_reasons.get(reason, 0) + 1
                continue

            # Use empty string if 'Abstract' is missing
            title_abstract = paper['ArticleTitle'] + " " + paper.get('Abstract', '')
            title_abstract_lower = title_abstract.lower()

            # Check for inclusion keywords
            if not any(keyword.lower() in title_abstract_lower for keyword in inclusion_keywords):
                reason = "Missing inclusion keywords"
                exclusion_reasons[reason] = exclusion_reasons.get(reason, 0) + 1
                continue

            # Check for exclusion keywords
            if any(keyword.lower() in title_abstract_lower for keyword in exclusion_keywords):
                reason = "Contains exclusion keywords"
                exclusion_reasons[reason] = exclusion_reasons.get(reason, 0) + 1
                continue

            # If all checks pass, include the paper
            filtered.append(paper)

        except Exception as e:
            logging.error(f"Error processing paper: {e}")

    if verbose:
        print("\nExclusion Reasons and Counts:")
        for reason, count in exclusion_reasons.items():
            print(f"- {reason}: {count}")

    return filtered

def save_results_to_csv(filtered_papers, file_name="filtered_papers.csv"):
    """
    Save filtered papers to a CSV file.
    
    Args:
        filtered_papers (list): A list of filtered paper dictionaries.
        file_name (str): The name of the CSV file to save.
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
        "(delusion OR persecution OR erotomanic OR cotard OR reference OR jealous OR capgras OR nihilistic OR somatic OR grandeur OR religious OR antichrist OR jerusalem) AND "
        "(brain OR MRI OR CT)"
    )
    inclusion_criteria = ["case report", "case series", "delusion", "psychosis", "brain", "mri", "ct"]
    exclusion_criteria = ["toxic", "drug-induced", "metabolic", "physiologic", "pharmacological"]

    # PubMed Search
    pubmed_results = search_pubmed(search_query, max_results=10000)  # Set to 10,000
    if pubmed_results:
        pubmed_articles = fetch_pubmed_details(pubmed_results, verbose=False)  # Ensure verbose is False
        print(f"\nNumber of PubMed articles retrieved: {len(pubmed_articles)}")
    else:
        pubmed_articles = []
        print("No PubMed articles found.")

    # Combine and Filter Articles
    articles = pubmed_articles
    print(f"Total number of articles retrieved: {len(articles)}")

    if articles:
        filtered_papers = filter_papers(
            articles,
            inclusion_keywords=inclusion_criteria,
            exclusion_keywords=exclusion_criteria,
            verbose=True  # Enable detailed exclusion messages
        )
        if filtered_papers:
            save_results_to_csv(filtered_papers)
            print(f"{len(filtered_papers)} papers match the criteria.")
        else:
            print("No papers matched the criteria.")
    else:
        print("No articles to process.")