__author__ = 'fbidu'

import datetime
import string
import Article
import ConfigParser
import sys


def has_special_char(text):
    for char in text:
        if char in string.punctuation:
            return True

    return False


def list_names(row, dates, columns):
    """Function that returns a dictionary with name:date for each row"""

    # List to store the skipped data
    skipped = []

    # Dictionary to store the names and the date they contacted us
    dates = dates

    # Looping through all the columns with author names
    for column in columns:

        name = string.capwords(row[column])

        try:
            name = name.decode('utf8')
        except UnicodeDecodeError:
            print "Decode error @ %s" % name

        # If there are punctuations or digits, we skip them
        if has_special_char(name) or len(name.split()) <= 1:
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


def get_articles_by_researcher(articles_ids, res, pubmedutils):
    """Function that gets all the articles in articles_ids"""

    # Looping through all the article IDs
    for article_id in articles_ids:
        # Getting the raw data from the article
        raw_article = pubmedutils.get_article(article_id)

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
                k = data['MedlineCitation']['KeywordList'][0]
            except (KeyError, IndexError):
                k = 0

            keywords = []

            if k:
                for keyword in k:
                    keywords.append(keyword)

            # Initializing an object
            article = Article.Article(pubmed_id=article_id,
                                      title=title,
                                      abstract=abstract,
                                      keywords=keywords)

            # Appending to the author listing
            res.articles.append(article)


def get_list_in_config(config_file, section_name, field_name, delimiter=','):
    values = []

    config = ConfigParser.ConfigParser()

    if len(config.read(config_file)) < 1:
        print "There are no settings on the supplied file."
        sys.exit(1)

    temp = config.get(section_name, field_name)

    if len(temp) < 1:
        return values
    else:
        values = [int(x) for x in temp.split(delimiter)]
        return values