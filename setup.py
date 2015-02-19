try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'medCrawl',
    'author': 'Felipe Bidu',
    'url': 'http://github.com/fbidu',
    'download_url': 'https://github.com/fbidu/medCrawl',
    'author_email': 'bidu.pub@gmail.com',
    'version': '0.1',
    'install_requires': ['nose', 'Biopython'],
    'packages': ['medcrawl'],
    'scripts': [],
    'name': 'medCrawl'
}

setup(**config)
__author__ = 'fbidu'
