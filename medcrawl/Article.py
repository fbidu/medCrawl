"""
This class provides the Article class
"""
__author__ = 'F. Bidu'


class Article(object):
    """This class holds the Article data"""
    pubmed_id = ""
    title = ""
    abstract = ""
    keywords = ""

    def __init__(self, pubmed_id, title, abstract, keywords):
        self.pubmed_id = str(pubmed_id)
        self.title = title
        self.abstract = abstract
        self.keywords = keywords

    def get_abstract(self):
        """Returns the article abstract"""
        return self.abstract

    def get_keywords(self):
        """Returns the article's keywords"""
        return self.keywords
