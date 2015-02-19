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

i = 0
confocal = []
dates = {}
print("Loading data")
with open('data/confocal.csv', 'rb') as csvfile:

    dataReader = csv.reader(csvfile, delimiter=',')
    print("CSV reading completed, parsing data")
    for row in dataReader:
        if i > 0:
            name = string.capwords(row[4])
            prof = string.capwords(row[1])
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

now = str(datetime.datetime.today().year)+'/'+str(datetime.datetime.today().month)

for item in dates:
    print item + ' ' + str(dates[item])

articles = []
Entrez.email = "Your.Name.Here@example.org"
for nome in confocal:
    print("Getting articles by "+nome)
    then = str(dates[nome].year)+'/'+str(dates[nome].month)
    handle = Entrez.esearch(db="pubmed", term=nome+"[AU] "+then+':'+now+'[DP]')
    record = Entrez.read(handle)
    handle.close()
    for id in record['IdList']:
        if id not in articles:
            articles.append(id)
print(str(len(articles))+" found!")
print(articles)
        # for data in getArticle(id):
        #     try:
        #         date = data['MedlineCitation']['Article']['ArticleDate'][0]['Day']+'/'+ \
        #            data['MedlineCitation']['Article']['ArticleDate'][0]['Month']+'/'+ \
        #            data['MedlineCitation']['Article']['ArticleDate'][0]['Year']
        #     except:
        #         pass
        #
        #     print('\t'+data['MedlineCitation']['Article']['ArticleTitle']+' Date: '+date)

# nm = open("names.lst", "w")
#
# for name in sorted_names:
#     nm.write(name[0]+': ')
#     nm.write(str(name[1]))
#     nm.write('\n')
#
# nm.close()

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
