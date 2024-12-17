# Paper Abstract Filterer

## Overview

**Paper Abstract Filterer** is a Python script designed to search, retrieve, filter, and save academic articles from PubMed based on specific inclusion and exclusion criteria. This tool leverages the [Biopython](https://biopython.org/) library to interact with the PubMed API, allowing researchers and professionals to efficiently gather and curate relevant literature for their studies or projects.

## Features

- **PubMed Search:** Perform comprehensive searches on PubMed using customized queries.
- **Batch Retrieval:** Fetch article details in batches of 100 to comply with NCBI rate limits.
- **Keyword Filtering:** 
  - **Inclusion Keywords:** Ensure articles contain specific keywords.
  - **Exclusion Keywords:** Exclude articles containing undesirable keywords.
- **Detailed Reporting:** Display reasons for exclusion and count of excluded articles.
- **CSV Export:** Save the filtered list of articles to a CSV file for easy analysis and sharing.
- **Progress Indicators:** Real-time updates on the fetching process to monitor progress without cluttering the console.

## Prerequisites

- **Python 3.6 or higher**
- **Internet Connection:** Required to access the PubMed API and download NLTK data.

## Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/yourusername/paper-abstract-filterer.git
   cd paper-abstract-filterer
   ```

2. **Create a Virtual Environment (Recommended):**

   ```bash
   python -m venv venv
   ```

3. **Activate the Virtual Environment:**

   - **Windows:**
     ```bash
     venv\Scripts\activate
     ```
   - **macOS/Linux:**
     ```bash
     source venv/bin/activate
     ```

4. **Install Required Dependencies:**

   ```bash
   pip install -r requirements.txt
   ```

   *If `requirements.txt` is not provided, install dependencies manually:*

   ```bash
   pip install biopython pandas nltk
   ```

5. **Download NLTK Data:**

   The script automatically checks and downloads the necessary NLTK data (`punkt` and `stopwords`) to the specified directory (`C:\Users\wooot\nltk_data`). Ensure that the script has the necessary permissions to create directories and download data.

## Configuration

1. **Set Your Email for Entrez:**

   PubMed requires users to provide an email address when using the Entrez API. Update the 

Entrez.email

 field in 

filter.py

 with your actual email address.

   ```python
   Entrez.email = "your.email@example.com"  # Replace with your actual email address
   ```

## Usage

1. **Define Your Search Query and Criteria:**

   The script uses predefined search queries, inclusion criteria, and exclusion criteria. Modify these variables in the `__main__` section of 

filter.py

 as needed.

   ```python
   search_query = (
       "(case reports[Publication Type] OR case series[Publication Type]) AND "
       "(delusion OR persecution OR erotomanic OR cotard OR reference OR jealous OR capgras OR nihilistic OR somatic OR grandeur OR religious OR antichrist OR jerusalem) AND "
       "(brain OR MRI OR CT)"
   )
   inclusion_criteria = ["case report", "case series", "delusion", "psychosis", "brain", "mri", "ct"]
   exclusion_criteria = ["toxic", "drug-induced", "metabolic", "physiologic", "pharmacological"]
   ```

2. **Run the Script:**

   Execute the script using Python:

   ```bash
   python filter.py
   ```

3. **Monitor Progress:**

   - The script will display real-time updates on fetching records, updating the status on a single line to keep the console output clean.
   - Upon completion, it will display the number of articles retrieved and the results of the filtering process.

4. **Output:**

   - **Console Output:**
     ```
     Punkt tokenizer found!
     Stopwords found!
     Searching PubMed for: (your search query)
     Number of articles found: 19177
     Note: Only the first 10000 articles will be fetched due to API limitations.
     Fetching records 1 to 100
     Fetching records 101 to 200
     ...
     Fetching records 9901 to 10000

     Number of PubMed articles retrieved: 10000
     Total number of articles retrieved: 10000

     Exclusion Reasons and Counts:
     - Missing required fields: 50
     - Missing inclusion keywords: 200
     - Contains exclusion keywords: 300

     Results saved to filtered_papers.csv
     9150 papers match the criteria.
     ```

   - **CSV File:**
     - The filtered results are saved to 

filtered_papers.csv

 in the same directory as the script. The CSV includes columns like `ArticleTitle`, `AuthorList`, `Abstract`, and `PubDate`.

## Troubleshooting

- **Syntax Errors:**
  - Ensure that all code blocks are properly indented. Python relies on indentation to define code structure.
  - Verify that all parentheses, brackets, and braces are correctly closed.

- **Function Argument Errors:**
  - If encountering errors like `TypeError: fetch_pubmed_details() got an unexpected keyword argument 'verbose'`, ensure that the function definition includes the `verbose` parameter.

- **NLTK Data Download Issues:**
  - Make sure the script has permission to create directories and download data to the specified `nltk_data_path`.
  - If manual intervention is needed, you can download NLTK data by running the following in a Python shell:

    ```python
    import nltk
    nltk.download('punkt', download_dir=r'C:\Users\wooot\nltk_data')
    nltk.download('stopwords', download_dir=r'C:\Users\wooot\nltk_data')
    ```

- **PubMed API Rate Limits:**
  - The script includes a 

time.sleep(0.4)

 between API calls to respect NCBI's rate limits. Avoid lowering this interval to prevent being blocked by the API.

## Customization

- **Adjust Batch Size and Rate Limiting:**
  - Modify the 

batch_size

 and 

time.sleep()

 duration in the `fetch_pubmed_details` function to better suit your needs or to comply with API policies.

- **Modify Inclusion and Exclusion Keywords:**
  - Update the `inclusion_criteria` and `exclusion_criteria` lists to tailor the filtering process to your specific requirements.

- **Change Output File Name:**
  - The default output file is 

filtered_papers.csv

. You can change this by providing a different `file_name` argument when calling `save_results_to_csv`.

## Dependencies

- [Biopython](https://biopython.org/)
- [Pandas](https://pandas.pydata.org/)
- [NLTK](https://www.nltk.org/)

Install dependencies using:

```bash
pip install biopython pandas nltk
```

## License

This project is licensed under the MIT License.

## Acknowledgements

- Thanks to the developers of Biopython, Pandas, and NLTK for providing robust libraries that make tasks like this possible.
- NCBI for providing the PubMed API, facilitating access to a vast repository of biomedical literature.

## Contact

For any questions or suggestions, please contact [James Thomas Rouse](mailto:james.thomas.rouse@gmail.com).