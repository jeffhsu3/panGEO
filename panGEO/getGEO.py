import urllib2, os, gzip, re
import pandas as pd

class geoQuery(object):
    base_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    pass

class BaseExpr(object):
    """ Base expression type
    """

    def __init__(self):
        pass


def _construct_geo_ftp(geo_id):
    """ Construct an ftp address from the GEO ID
    """

    url_header = 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/'
    geo_type = geo_id[0:3]
    if geo_type == 'GSE':
        url_header += 'SOFT/by_series/' + geo_id +\
                '/%s_family.soft.gz' % (geo_id)
    elif geo_type == 'GPL':
        #:TODO finish this
        pass
    return(url_header)

def _parse_meta_data():
    pass


def getGEO(geo_id, out_file='/tmp/'):
    temp = urllib2.urlopen(_construct_geo_ftp(geo_id))
    # :TODO this needs to not be OS specific
    # :TODO seems terribly inefficient to save and reload the gzip.
    # see if you can refactor to use IOStream.

    all_data = {}
    out = os.path.join("/tmp", geo_id)
    open(out, 'wb').write(temp.read())
    temp = gzip.open(out, 'rb')
    content = temp.read()
    temp.close()
    parseGSE(content)



def loadGEO(GEO_file):
    """ Loads a locally downloaded GEO file
    """
    temp = gzip.open(GEO_file, 'rb')
    content = temp.read()
    return(content)


def parseGSE(content):
    """ All in memory based right now how can we change?  Liberal use of
    seeks?
    """

    gse = {}

    gse_platforms = re.compile(u'\\^PLATFORM')
    gse_samples = re.compile(u'\\^SAMPLE')


    platforms = gse_platforms.finditer(content)
    samples = gse_samples.finditer(content)

    gsm ={}
    gsm['starts'] = []
    gsm['gsm_id'] = []

    samplestarts = []
    for i in samples:
        endline = content.find('\n', i.start())
        gsm['starts'].append(i.start())
        gsm['gsm_id'].append(content[i.start():endline].split(' = ')[1])
    print(gsm['starts'])
    chr_to_read = [gsm['starts'][i+1]-1 for i\
                         in range(len(gsm['starts'])-1)]
    # Assuming the last one ends at the end of the file, also -1 for EOF
    chr_to_read.append(len(content)-1)
    print(len(gsm['starts']))

    for s, e in zip(gsm['starts'], chr_to_read):
        k = content[s:e]
        re_plat_id = re.compile(u'\\!Sample_platform_id.*$', re.MULTILINE)
        plat_id = re_plat_id.search(k).group().split(" = ")[1]
        re_data_start = re.compile(u'_table_begin$', re.MULTILINE)
        data_start = re_data_start.search(k).end() + 1
        data = k[data_start:].split("\n")
        print(data[-1])

class eSet(object):
    """ A python version of bioconductor's eset using panda's dataframes
    """
    pass
