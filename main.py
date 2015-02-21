"""
Script that crawls Pubmed for articles
by some set of authors after some date
"""

__author__ = 'F. Bidu'

import sys
import csv
import datetime
import string
import glob

from medcrawl import Article, Researcher, PubmedUtils

from jinja2 import Environment, FileSystemLoader

# Initializing the pubmed utils
PU = PubmedUtils.PubMedUtils(email="bidu.pub@gmail.com")

# Possible modes: namesOnly, idsOnly, fullArticles
MODE = "fullArticles"

# debug mode (2), normal mode (1), silent mode (0)
DEBUG = 2

# List with the columns with data to be gathered
COLUMN_WITH_AUTHORS = [1, 4]

# List with the csv files at the data folder
FILES = glob.glob("data/*.csv")

# If there are no files in the data folder, the program exits
if len(FILES) < 1:
    print "No CSV files found at data folder"
    sys.exit(1)

# Setting the templates folder
TEMPLATE_FOLDER = '.'

# Loading the template environment
ENV = Environment(loader=FileSystemLoader(TEMPLATE_FOLDER))


def list_names(row, dates):
    """Function that returns a dictionary with name:date for each row"""

    # List to store the skipped data
    skipped = []

    # Dictionary to store the names and the date they contacted us
    dates = dates

    # Looping through all the columns with author names
    for column in COLUMN_WITH_AUTHORS:

        name = string.capwords(row[column])

        # If there are punctuations or digits, we skip them
        if (string.punctuation in name) or \
                (string.digits in name) or \
                (len(name.split()) <= 1):
            skipped.append(name)
        else:
            # If the name is valid,
            # we'll try to find the first date it occurred
            if row[0] == '':
                continue

            current_date = datetime.datetime.strptime(
                row[0],
                "%d/%m/%Y %H:%M:%S").date()

            # If the name is already on the dictionary...
            if name in dates:
                # We just check for the date
                if dates[name] > current_date:
                    dates[name] = current_date
            else:
                dates[name] = current_date

    return dates


def get_articles_by_researcher(articles_ids, res):
    """Function that gets all the articles in articles_ids"""

    # Looping through all the article IDs
    for article_id in articles_ids:
        # Getting the raw data from the article
        raw_article = PU.get_article(article_id)

        # Looping through the acquired data
        for data in raw_article:
            # Extracting the MedLine[Article] section
            md_article = data['MedlineCitation']['Article']

            # Extracting the title
            title = md_article['ArticleTitle']

            try:
                title = title.decode('utf8')
            except UnicodeEncodeError:
                pass

            # Extracting the abstract
            try:
                abstract = md_article['Abstract']['AbstractText'][0]
                try:
                    abstract = abstract.decode('utf8')
                except UnicodeEncodeError:
                    pass

            except KeyError:
                abstract = "Abstract not available"

            # Extracting the keywords
            try:
                keywords = 'oi, oi, oi'
                # keywords = data['MedlineCitation']['KeywordList']
            except KeyError:
                keywords = "Keywords not available"

            # Initializing an object
            article = Article.Article(pubmed_id=article_id,
                                      title=title,
                                      abstract=abstract,
                                      keywords=keywords)

            # Appending to the author listing
            res.articles.append(article)


def main():
    """
    Main function
    """

    print "Loading data..."

    # Looping through the list of files
    for csv_file in FILES:

        # Opening the file
        print "Opening the file %s" % csv_file

        with open("%s" % csv_file, 'rb') as data_file:
            # Feeding to the reader
            data_reader = csv.reader(data_file, delimiter=',')

            # Skipping the first row - contains only column names
            data_reader.next()

            print "File %s successfully read! Parsing data" % csv_file

            dates = {}

            # Looping through all the rows of the file
            for row in data_reader:
                dates.update(list_names(row, dates))

            names = dates.keys()
            names.sort()

            # List to hold all of the researchers
            researchers = []

            # Looping through all the names
            for researcher in names:

                articles_ids = []

                # Instantiating a new researcher
                res = Researcher.Researcher(name=researcher)

                print "Getting articles by " + researcher if DEBUG else 0

                # Locating articles by date
                then = "%s/%s" % (str(dates[researcher].year),
                                  str(dates[researcher].month))

                articles_ids += PU.get_articles_in_range_by_author(researcher,
                                                              range_start=then)

                # getting data from all of the articles
                # if the mode is 'fullArticles'
                if MODE == 'fullArticles':
                    get_articles_by_researcher(articles_ids, res)

                if len(res.articles) > 0:
                    researchers.append(res)

            # Loading the template
            template = ENV.get_template('template.html')

            # Generating the output
            output = template.render(researchers=researchers).encode("utf8")

            # Saving the results
            with open("%s.html" % csv_file[:-4], "wb") as output_file:
                output_file.write(output)
    sys.exit(0)

if __name__ == '__main__':
    main()
