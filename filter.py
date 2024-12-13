from Bio import Entrez
import pandas as pd
import nltk
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
import string
import os
import time  # For rate limiting

# Define the custom nltk_data path
nltk_data_path = r'C:\Users\wooot\nltk_data'
os.environ['NLTK_DATA'] = nltk_data_path
os.makedirs(nltk_data_path, exist_ok=True)
nltk.data.path.append(nltk_data_path)

# Debug: Print the NLTK search paths
print("NLTK search paths:", nltk.data.path)

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

def search_pubmed(query, max_results=1000):
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

def fetch_pubmed_details(search_results):
    """
    Fetch detailed PubMed information using WebEnv and QueryKey.
    """
    records = []
    try:
        count = int(search_results["Count"])
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]
        batch_size = 100
        for start in range(0, count, batch_size):
            end = min(start + batch_size, count)
            print(f"Fetching records {start+1} to {end}")
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
                records.extend(result['PubmedArticle'])
            else:
                print("No 'PubmedArticle' found in the result.")
                print("Result keys:", result.keys())
        return records
    except Exception as e:
        print(f"An error occurred while fetching PubMed details: {e}")
        return records

def filter_papers(papers, inclusion_criteria, exclusion_criteria):
    """
    Filter papers based on inclusion and exclusion criteria.
    """
    filtered = []
    for paper in papers:
        # Type check to ensure 'paper' is a dictionary
        if not isinstance(paper, dict):
            print(f"Skipping paper; unexpected type: {type(paper)}")
            continue
        try:
            # Extract title and abstract
            medline_citation = paper.get("MedlineCitation", {})
            if not medline_citation:
                print("Skipping paper; 'MedlineCitation' not found")
                continue
            article = medline_citation.get("Article", {})
            title = article.get("ArticleTitle", "")
            abstract_sections = article.get("Abstract", {}).get("AbstractText", [])
            
            # Combine abstract sections if necessary
            if isinstance(abstract_sections, list):
                abstract = " ".join(abstract_sections)
            elif isinstance(abstract_sections, str):
                abstract = abstract_sections
            else:
                abstract = ""
            
            full_text = (title + " " + abstract).lower()
            
            # Apply inclusion criteria
            has_pub_type = any(pt in full_text for pt in ["case report", "case series"])
            has_delusion = "delusion" in full_text
            has_imaging = any(im in full_text for im in ["brain scan", "mri", "ct"])
            is_included = has_pub_type and has_delusion and has_imaging

            # Apply exclusion criteria
            is_excluded = any(conf.lower() in full_text for conf in exclusion_criteria)

            if is_included and not is_excluded:
                pmid_element = medline_citation.get("PMID", "")
                # Access the PMID value directly
                if hasattr(pmid_element, '_value'):
                    pmid = pmid_element._value
                else:
                    pmid = str(pmid_element)
                
                filtered.append({
                    "Title": title,
                    "Abstract": abstract,
                    "ID": pmid,
                })
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
    # Define your search term and criteria
    search_query = (
        "(case reports[Publication Type] OR case series[Publication Type]) AND "
        "(delusion) AND "
        "(brain scan OR MRI OR CT)"
    )
    inclusion_criteria = ["case report", "case series", "delusion", "brain scan", "mri", "ct"]
    exclusion_criteria = ["toxic", "drug induced", "metabolic", "physiologic", "pharmacological"]

    # Perform search and fetch details
    search_results = search_pubmed(search_query, max_results=1000)
    if not search_results:
        print("No results found.")
    else:
        id_list = search_results["IdList"]
        print(f"Number of article IDs retrieved: {len(id_list)}")
        articles = fetch_pubmed_details(search_results)
        if not articles:
            print("No articles found.")
        else:
            print(f"Number of articles retrieved: {len(articles)}")
            # Filter papers
            filtered_papers = filter_papers(articles, inclusion_criteria, exclusion_criteria)
            if not filtered_papers:
                print("No papers matched the criteria.")
            else:
                # Save results
                save_results_to_csv(filtered_papers)
                print(f"{len(filtered_papers)} papers match the criteria.")