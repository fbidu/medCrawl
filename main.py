__author__ = 'fbidu'

import sys
import csv
import time

try:
    from Bio import Entrez
except ImportError:
    print("Biopython not found!")
    sys.exit(1)

print("Loading data")

def getArticle(pubmedid):
    fetcher = Entrez.efetch(db="pubmed", id=pubmedid, retmode="xml")
    return Entrez.parse(fetcher)

i = 0
confocal = []
with open('data/confocal.csv', 'rb') as csvfile:

    dataReader = csv.reader(csvfile, delimiter=',')
    for row in dataReader:
        if i > 0:
            if row[0] != '':
                print(
                    "year = " + str(time.strptime(row[0], "%d/%m/%Y %H:%M:%S").tm_year) +
                    " month = " + str(time.strptime(row[0], "%d/%m/%Y %H:%M:%S").tm_mon) +
                    " day = " + str(time.strptime(row[0], "%d/%m/%Y %H:%M:%S").tm_mday) +
                    " email = " + row[4]
                )
        i += 1

Entrez.email = "Your.Name.Here@example.org"
handle = Entrez.esearch(db="pubmed", retmax=10, term="Sandra Krauchenco[AU]")
record = Entrez.read(handle)
handle.close()

for id in record['IdList']:
    for data in getArticle(id):
        print(data['MedlineCitation']['Article']['ArticleTitle'])
