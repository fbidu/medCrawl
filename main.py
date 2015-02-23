"""
Script that crawls Pubmed for articles
by some set of authors after some date
"""

__author__ = 'F. Bidu'

import sys
import csv
import glob
import ConfigParser
import os

from medcrawl import Researcher, PubmedUtils, Utils

from jinja2 import Environment, FileSystemLoader

# Initializing the pubmed utils
PU = PubmedUtils.PubMedUtils(email="bidu.pub@gmail.com")

# Loading the settings file
CONF = ConfigParser.ConfigParser()

CONFIG_FILE = 'medcrawl.conf'

if len(CONF.read(CONFIG_FILE)) < 1:
    print "No config file found!"
    sys.exit(1)

# Possible modes: namesOnly, idsOnly, fullArticles
MODE = CONF.get('medCrawl', 'mode')

# debug mode (2), normal mode (1), silent mode (0)
DEBUG = CONF.getint('medCrawl', 'debug')

# Getting the data folder
DATA_FOLDER = CONF.get('medCrawl', 'dataFolder')

# Getting the output folder
OUTPUT_FOLDER = CONF.get('medCrawl', 'outputFolder')

# List with the csv files at the data folder
FILES = glob.glob("%s/*.csv" % DATA_FOLDER)

# Getting the templates folder
TEMPLATE_FOLDER = CONF.get('medCrawl', 'templateFolder')

# If there are no files in the data folder, the program exits
if len(FILES) < 1:
    print "No CSV files found at data folder"
    sys.exit(1)

# Loading the template environment
ENV = Environment(loader=FileSystemLoader(TEMPLATE_FOLDER))


def main():
    """
    Main function
    """

    print "Loading data..."

    # Looping through the list of files
    for csv_file in FILES:

        # Getting the base file name
        base_file_name = os.path.basename(csv_file)[:-4]

        # List with the columns with data to be gathered
        column_with_authors = Utils.get_list_in_config(CONFIG_FILE, base_file_name, 'nameColumns')

        # Opening the file
        print "Opening the file %s" % csv_file

        with open("%s" % csv_file, 'r') as data_file:
            # Feeding to the reader
            data_reader = csv.reader(data_file, delimiter=',')

            # Skipping the first row - contains only column names
            data_reader.next()

            print "File %s successfully read! Parsing data" % csv_file

            dates = {}

            # Looping through all the rows of the file
            for row in data_reader:
                dates.update(Utils.list_names(row, dates, column_with_authors))

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
                    Utils.get_articles_by_researcher(articles_ids, res, PU)

                if len(res.articles) > 0:
                    researchers.append(res)

            # Loading the template
            template = ENV.get_template('template.html')

            # Generating the output
            output = template.render(researchers=researchers, title=base_file_name).encode("utf8")

            # Saving the results
            with open("%s/%s.html" % (OUTPUT_FOLDER, base_file_name), "wb") as output_file:
                output_file.write(output)
    sys.exit(0)


if __name__ == '__main__':
    main()
