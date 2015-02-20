__author__ = 'fbidu'

import sys
import csv
import datetime
import string
import glob
from jinja2 import Environment, FileSystemLoader

try:
    from Bio import Entrez
except ImportError:
    print("Biopython not found!")
    sys.exit(1)


class Article:
    id = ""
    title = ""
    abstract = ""
    keywords = ""

    def __init__(self, id, title, abstract, keywords):
        self.id = str(id)
        self.title = title
        self.abstract = abstract
        self.keywords = keywords


class Researcher:
    name = ""
    articles = []
    articlesByKeywords = {}

    def __init__(self, name, articles = []):
        self.name = name
        self.articles = articles


def getArticle(pubmedid):
    fetcher = Entrez.efetch(db="pubmed", id=pubmedid, retmode="xml")
    return Entrez.parse(fetcher)


def getArticlesInRangeByAuthor(author, rangeStart, rangeEnd='', terms=''):
    # if there's no rangeEnd set, we'll set it to the current year/month
    if rangeEnd == '':
        rangeEnd = str(datetime.datetime.today().year) + '/' + str(datetime.datetime.today().month)

    Entrez.email = "Your.Name.Here@example.org"
    term = "%s [AU] %s:%s[DP] %s" % (author, rangeStart, rangeEnd, terms)
    print(term)
    handle = Entrez.esearch(db="pubmed", term=term)
    record = Entrez.read(handle)
    handle.close()

    articles = []
    for articleId in record['IdList']:
        if articleId not in articles:
            articles.append(articleId)

    return articles


# Possible modes: namesOnly, idsOnly, fullArticles
mode = "fullArticles"

# debug mode (2), normal mode (1), silent mode (0)
debug = 1


# List with the columns with data to be gathered
columnWithAuthors = [1, 4]

# List with the csv files at the data folder
files = glob.glob("data/*.csv")

# If there are no files in the data folder, the program exits
if len(files) < 1:
    print("No CSV files found at data folder")
    sys.exit(1)

print("Loading data...")

# Looping through the list of files
for csvFile in files:

    # Opening the file
    with open("%s" % csvFile, 'rb') as dataFile:
        # Feeding to the reader
        dataReader = csv.reader(dataFile, delimiter=',')

        # Skipping the first row - contains only column names
        dataReader.next()

        # List to store the skipped data
        skipped = []

        # Dictionary to store the names and the date they contacted us
        dates = {}

        print("File %s sucessfully read! Parsing data" % csvFile)

        # Looping through all the rows of the file
        for row in dataReader:

            # Looping through all the columns with author names
            for column in columnWithAuthors:

                name = string.capwords(row[column])

                # If there are punctuations or digits, we skip them
                if (string.punctuation in name) or (string.digits in name) or (len(name.split()) <= 1):
                    skipped.append(name)
                else:
                    # If the name is valid, we'll try to find the first date it occurred
                    currentDate = datetime.datetime.strptime(row[0], "%d/%m/%Y %H:%M:%S").date()

                    # If the name is already on the dictionary...
                    if name in dates:
                        # We just check for the date
                        if dates[name] > currentDate:
                            dates[name] = currentDate
                    else:
                        dates[name] = currentDate

        # We've now read the entire file, so we'll start getting the data
        names = dates.keys()
        names.sort()

        # List to hold all of the researchers
        researchers = []

        # Looping through all the names
        for researcher in names:

            articlesIds = []
            # Instantiating a new researcher
            res = Researcher(name=researcher.decode("utf8"))

            print("Getting articles by " + researcher) if debug else 0

            # Locating articles by date
            then = str(dates[researcher].year) + '/' + str(dates[researcher].month)
            articlesIds += getArticlesInRangeByAuthor(researcher, rangeStart=then)

            # getting data from all of the articles if the mode is 'fullArticles'
            if mode == 'fullArticles':

                # Looping through all the article IDs
                for articleId in articlesIds:
                    # Getting the raw data from the article
                    rawArticle = getArticle(articleId)

                    for data in rawArticle:
                        # Initiating on an object

                        title = data['MedlineCitation']['Article']['ArticleTitle'].decode('utf8')

                        try:
                            abstract = data['MedlineCitation']['Article']['Abstract']['AbstractText'][0].decode("utf8")
                        except:
                            abstract = "Abstract not available"

                        try:
                            keywords = 'oi, oi, oi'
                            # keywords = data['MedlineCitation']['KeywordList']
                        except:
                            keywords = "Keywords not available"

                        article = Article(id=articleId,
                                          title=title,
                                          abstract=abstract,
                                          keywords=keywords)

                        # Appending to the author listing
                        res.articles.append(article)

            researchers.append(res)

        env = Environment(loader=FileSystemLoader('.'))
        template = env.get_template('template.html')
        output_from_parsed_template = template.render(researchers=researchers).encode("utf8")
        # print output_from_parsed_template

        # to save the results
        with open("my_new_file.html", "wb") as fh:
            fh.write(output_from_parsed_template)

        sys.exit(0)

# print("Saving list in file")
# nm = open("confocalView.lst", "w")
# nm.writelines("""<html>
# 	<head>
# 		<title>Articles that might be related to Confocal</title>
# 		<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.2/css/bootstrap.min.css">
# 	</head>
# 	<body>
# 		<div class="jumbotron">
#       <div class="container">""")
#
#
# articles += getArticlesInRangeByAuthor(name, rangeStart=then)
# microscopy += getArticlesInRangeByAuthor(name, rangeStart=then, terms="microscopy[ALL]")
# confocalk += getArticlesInRangeByAuthor(name, rangeStart=then, terms="confocal[ALL]")
# if len(articles) > 0:
#     nm.write("<h2>Articles by %s</h2>" % name.decode('utf8'))
#
#     if len(microscopy) > 0:
#         nm.write("<h3>Articles containing 'microscopy'</h3>")
#         for article in microscopy:
#
#             nm.writelines("""<div class="panel panel-danger"> <div class="panel-heading">""")
#             for data in getArticle(article):
#                 nm.write("<a href=http://www.ncbi.nlm.nih.gov/pubmed/?term=%s>%s</a>" % (
#                     str(article), data['MedlineCitation']['Article']['ArticleTitle'].encode('utf8')))
#                 nm.write("</div>")
#                 nm.write("<div class='panel-body'>")
#                 try:
#                     nm.write("<small>%s</small>" % str(
#                         data['MedlineCitation']['Article']['Abstract']['AbstractText'][0].encode("utf8")))
#                 except:
#                     nm.write("<small>No abstract available!</small>")
#
#                 nm.write("</div>")
#                 nm.write('<div class="panel-footer">')
#                 for keys in data['MedlineCitation']['KeywordList']:
#                     for keyword in keys:
#                         nm.write(keyword.encode('utf-8') + '; ')
#                 nm.write('</div></div>')
#     articles = list(set(articles) - set(microscopy))
#
#     if len(confocalk) > 0:
#         nm.write("<h3>Articles containing 'confocal'</h3>")
#         for article in confocalk:
#
#             nm.writelines("""<div class="panel panel-danger"> <div class="panel-heading">""")
#             for data in getArticle(article):
#                 nm.write("<a href=http://www.ncbi.nlm.nih.gov/pubmed/?term=%s>%s</a>" % (
#                     str(article), data['MedlineCitation']['Article']['ArticleTitle'].encode('utf8')))
#                 nm.write("</div>")
#                 nm.write("<div class='panel-body'>")
#                 try:
#                     nm.write("<small>%s</small>" % str(
#                         data['MedlineCitation']['Article']['Abstract']['AbstractText'][0].encode("utf8")))
#                 except:
#                     nm.write("<small>No abstract available!</small>")
#
#                 nm.write("</div>")
#                 nm.write('<div class="panel-footer">')
#                 for keys in data['MedlineCitation']['KeywordList']:
#                     for keyword in keys:
#                         nm.write(keyword.encode('utf-8') + '; ')
#                 nm.write('</div></div>')
#     articles = list(set(articles) - set(confocal))
#
#     for article in articles:
#         nm.writelines("""<div class="panel panel-default"> <div class="panel-heading">""")
#         for data in getArticle(article):
#             nm.write("<a href=http://www.ncbi.nlm.nih.gov/pubmed/?term=%s>%s</a>" % (
#                 str(article), data['MedlineCitation']['Article']['ArticleTitle'].encode('utf8')))
#             nm.write("</div>")
#             nm.write("<div class='panel-body'>")
#             try:
#                 nm.write("<small>%s</small>" % str(
#                     data['MedlineCitation']['Article']['Abstract']['AbstractText'][0].encode("utf8")))
#             except:
#                 nm.write("<small>No abstract available!</small>")
#
#             nm.write("</div>")
#             nm.write('<div class="panel-footer">')
#             for keys in data['MedlineCitation']['KeywordList']:
#                 for keyword in keys:
#                     nm.write(keyword.encode('utf-8') + '; ')
#             nm.write('</div></div>')
#
# nm.close()
# # for data in getArticle(id):
# # try:
# # date = data['MedlineCitation']['Article']['ArticleDate'][0]['Day']+'/'+ \
# # data['MedlineCitation']['Article']['ArticleDate'][0]['Month']+'/'+ \
# # data['MedlineCitation']['Article']['ArticleDate'][0]['Year']
# # except:
# # pass
# #
# # print('\t'+data['MedlineCitation']['Article']['ArticleTitle']+' Date: '+date)
#
#
# # if row[0] != '':
# # print(
# #         "year = " + str(time.strptime(row[0], "%d/%m/%Y %H:%M:%S").tm_year) +
# #         " month = " + str(time.strptime(row[0], "%d/%m/%Y %H:%M:%S").tm_mon) +
# #         " day = " + str(time.strptime(row[0], "%d/%m/%Y %H:%M:%S").tm_mday) +
# #         " email = " + row[4]
# #     )
# # Entrez.email = "Your.Name.Here@example.org"
# # handle = Entrez.esearch(db="pubmed", retmax=10, term="Sandra Krauchenco[AU]")
# # record = Entrez.read(handle)
# # handle.close()
# #
# # for id in record['IdList']:
# #     for data in getArticle(id):
# #         print(data['MedlineCitation']['Article']['ArticleTitle'])
