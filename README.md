# Paper Abstract Filterer

This Python script searches PubMed for scientific articles based on a user-defined query and processes the retrieved abstracts for further analysis. It leverages Biopython's Entrez module to interact with the PubMed database and NLTK for natural language processing tasks.

## Features

- **PubMed Search**: Connects to PubMed and searches for articles matching the specified query.
- **Abstract Retrieval**: Fetches abstracts of the found articles for analysis.
- **Text Processing**: Utilizes NLTK to tokenize text and remove stopwords, preparing the data for downstream tasks.

## Requirements

- Python 3.x
- [Biopython](https://biopython.org/)
- [NLTK](https://www.nltk.org/)

## Setup

1. **Install Dependencies**:

   ```bash
   pip install biopython nltk
   ```

2. **Download NLTK Data Files**:

   The script checks for the necessary NLTK data files (`punkt` and 

stopwords

) and downloads them if they are not present.

## Usage

1. **Set Your Email**:

   Modify the 

Entrez.email

 variable in 

filter.py

 to include your email address. This is required by NCBI to identify users.

2. **Run the Script**:

   ```bash
   python filter.py
   ```

3. **Enter Your Query**:

   When prompted, input your PubMed search query.

## Notes

- The script includes rate limiting to comply with NCBI's usage policies.
- NLTK data files are stored in a custom directory specified by 

nltk_data_path

.

## License

This project is licensed under the MIT License.