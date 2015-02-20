"""
This module provides only the Researcher class
"""
__author__ = 'F. Bidu'


class Researcher(object):
    """
    This class holds the data of a researcher
    """
    name = ""
    articles = []
    articlesByKeywords = {}

    def __init__(self, name, articles=None):
        articles = articles or []
        self.name = name
        self.articles = articles

    def list_articles(self):
        """Lists all the articles"""
        return self.articles

    def list_articles_by_keyword(self, keyword):
        """Lists all the articles within a keyword"""
        return self.articlesByKeywords[keyword]
