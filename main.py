__author__ = 'fbidu'

import sys
import csv
import datetime
import string

try:
    from Bio import Entrez
except ImportError:
    print("Biopython not found!")
    sys.exit(1)


def dumpclean(obj):
    if type(obj) == dict:
        for k, v in obj.items():
            if hasattr(v, '__iter__'):
                print k
                dumpclean(v)
            else:
                print '%s : %s' % (k, v)
    elif type(obj) == list:
        for v in obj:
            if hasattr(v, '__iter__'):
                dumpclean(v)
            else:
                print v
    else:
        print obj


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


i = 0
confocal = []
dates = {}
print("Loading data")
with open('data/confocal.csv', 'rb') as csvfile:
    dataReader = csv.reader(csvfile, delimiter=',')
    print("CSV reading completed, parsing data")
    for row in dataReader:
        if i > 0:
            name = string.capwords(row[4]).decode("utf8")
            prof = string.capwords(row[1]).decode("utf8")
            if '/' in name or name == '':
                pass
            else:
                date = datetime.datetime.strptime(row[0], "%d/%m/%Y %H:%M:%S").date()

                if name in dates:
                    if dates[name] > date:
                        dates[name] = date
                else:
                    dates[name] = date
                    if len(name.split(' ')) > 1:
                        confocal.append(name)

                if prof in dates:
                    if dates[prof] > date:
                        dates[prof] = date
                else:
                    dates[prof] = date
                    if len(prof.split(' ')) > 1:
                        confocal.append(prof)
        i += 1

sorted_names = sorted(dates.items())

now = str(datetime.datetime.today().year) + '/' + str(datetime.datetime.today().month)

print("Saving list in file")
nm = open("confocalView.lst", "w")
nm.writelines("""<html>
	<head>
		<title>Articles that might be related to Confocal</title>
		<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.2/css/bootstrap.min.css">
	</head>
	<body>
		<div class="jumbotron">
      <div class="container">""")

confocal.sort()

for name in confocal:
    articles = []
    microscopy = []
    confocalk = []
    print("Getting articles by " + name)
    then = str(dates[name].year) + '/' + str(dates[name].month)
    articles += getArticlesInRangeByAuthor(name, rangeStart=then)
    microscopy += getArticlesInRangeByAuthor(name, rangeStart=then, terms="microscopy[ALL]")
    confocalk += getArticlesInRangeByAuthor(name, rangeStart=then, terms="confocal[ALL]")
    if len(articles) > 0:
        nm.write("<h2>Articles by %s</h2>" % name.decode('utf8'))

        if len(microscopy) > 0:
            nm.write("<h3>Articles containing 'microscopy'</h3>")
            for article in microscopy:

                nm.writelines("""<div class="panel panel-danger"> <div class="panel-heading">""")
                for data in getArticle(article):
                    nm.write("<a href=http://www.ncbi.nlm.nih.gov/pubmed/?term=%s>%s</a>" % (
                        str(article), data['MedlineCitation']['Article']['ArticleTitle'].encode('utf8')))
                    nm.write("</div>")
                    nm.write("<div class='panel-body'>")
                    try:
                        nm.write("<small>%s</small>" % str(
                            data['MedlineCitation']['Article']['Abstract']['AbstractText'][0].encode("utf8")))
                    except:
                        nm.write("<small>No abstract available!</small>")

                    nm.write("</div>")
                    nm.write('<div class="panel-footer">')
                    for keys in data['MedlineCitation']['KeywordList']:
                        for keyword in keys:
                            nm.write(keyword.encode('utf-8') + '; ')
                    nm.write('</div></div>')
        articles = list(set(articles) - set(microscopy))

        if len(confocalk) > 0:
            nm.write("<h3>Articles containing 'confocal'</h3>")
            for article in confocalk:

                nm.writelines("""<div class="panel panel-danger"> <div class="panel-heading">""")
                for data in getArticle(article):
                    nm.write("<a href=http://www.ncbi.nlm.nih.gov/pubmed/?term=%s>%s</a>" % (
                        str(article), data['MedlineCitation']['Article']['ArticleTitle'].encode('utf8')))
                    nm.write("</div>")
                    nm.write("<div class='panel-body'>")
                    try:
                        nm.write("<small>%s</small>" % str(
                            data['MedlineCitation']['Article']['Abstract']['AbstractText'][0].encode("utf8")))
                    except:
                        nm.write("<small>No abstract available!</small>")

                    nm.write("</div>")
                    nm.write('<div class="panel-footer">')
                    for keys in data['MedlineCitation']['KeywordList']:
                        for keyword in keys:
                            nm.write(keyword.encode('utf-8') + '; ')
                    nm.write('</div></div>')
        articles = list(set(articles) - set(confocal))

        for article in articles:
            nm.writelines("""<div class="panel panel-default"> <div class="panel-heading">""")
            for data in getArticle(article):
                nm.write("<a href=http://www.ncbi.nlm.nih.gov/pubmed/?term=%s>%s</a>" % (
                    str(article), data['MedlineCitation']['Article']['ArticleTitle'].encode('utf8')))
                nm.write("</div>")
                nm.write("<div class='panel-body'>")
                try:
                    nm.write("<small>%s</small>" % str(
                        data['MedlineCitation']['Article']['Abstract']['AbstractText'][0].encode("utf8")))
                except:
                    nm.write("<small>No abstract available!</small>")

                nm.write("</div>")
                nm.write('<div class="panel-footer">')
                for keys in data['MedlineCitation']['KeywordList']:
                    for keyword in keys:
                        nm.write(keyword.encode('utf-8') + '; ')
                nm.write('</div></div>')

nm.close()
# for data in getArticle(id):
# try:
# date = data['MedlineCitation']['Article']['ArticleDate'][0]['Day']+'/'+ \
# data['MedlineCitation']['Article']['ArticleDate'][0]['Month']+'/'+ \
# data['MedlineCitation']['Article']['ArticleDate'][0]['Year']
# except:
# pass
#
#     print('\t'+data['MedlineCitation']['Article']['ArticleTitle']+' Date: '+date)


# if row[0] != '':
#     print(
#         "year = " + str(time.strptime(row[0], "%d/%m/%Y %H:%M:%S").tm_year) +
#         " month = " + str(time.strptime(row[0], "%d/%m/%Y %H:%M:%S").tm_mon) +
#         " day = " + str(time.strptime(row[0], "%d/%m/%Y %H:%M:%S").tm_mday) +
#         " email = " + row[4]
#     )
# Entrez.email = "Your.Name.Here@example.org"
# handle = Entrez.esearch(db="pubmed", retmax=10, term="Sandra Krauchenco[AU]")
# record = Entrez.read(handle)
# handle.close()
#
# for id in record['IdList']:
#     for data in getArticle(id):
#         print(data['MedlineCitation']['Article']['ArticleTitle'])
