"""
This module provides de PubmedUtils class
"""
__author__ = 'F. Bidu'

import sys
from datetime import datetime

try:
    from Bio import Entrez as E
except ImportError:
    print "Biopython not found!"
    sys.exit(1)


class PubMedUtils(object):
    """
    This class provides simple functions to manipulate the
    Pubmed API via Entrez
    """
    email = "Your.Name.Here@example.org"

    def __init__(self, email):
        self.email = email

    @classmethod
    def get_article(cls, pubmed_id):
        """
        This function returns the full data of the
        article id'd by pubmed_id
        """

        fetcher = E.efetch(db="pubmed", id=pubmed_id, retmode="xml")
        return E.parse(fetcher)

    def get_articles_in_range_by_author(self, author,
                                        range_start, range_end='',
                                        terms=''):
        """
        This function returns all of the article IDs
        of the articles a given 'author' published
        between 'range_start' and 'range_end'
        Optionally, additional terms may be provided
        """
        # if there's no rangeEnd set, we'll set it to the current year/month
        if range_end == '':
            range_end = "%s/%s" % (str(datetime.today().year),
                                   str(datetime.today().month))

        E.email = self.email
        term = "%s [AU] %s:%s[DP] %s" % (author, range_start, range_end, terms)

        handle = E.esearch(db="pubmed", term=term)
        record = E.read(handle)
        handle.close()

        articles = []
        for article_id in record['IdList']:
            if article_id not in articles:
                articles.append(article_id)
        print "got %d articles" % len(articles)
        return articles
