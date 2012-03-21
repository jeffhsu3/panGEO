import urllib2, os, gzip, re, tempfile, cStringIO, StringIO, collections
import pandas as pd


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


def parseGSM(fileobj, size, chunksize = 100000):
    seen = 0
    first_line = fileobj.readline()
    seen += len(first_line)
    acession = first_line.strip("\n").split(" = ")[1]
    meta_data = {}
    # This is slow
    while True:
        line = fileobj.readline()
        seen += len(line)
        if re.search(u'\\!Sample_title', line):
            meta_data['sample'] = line.rstrip('\n').split(" = ")[1]
        elif '!Sample_platform_id' in line:
            meta_data['platform'] = line.rstrip('\n').split(" = ")[1]
        else:
            pass
        if re.search('_table_begin', line):
            start = fileobj.tell()
            break
        else: pass
    # Get the line counts to read
    if size:
        data = fileobj.read(size-seen)
    else:
        data = fileobj.read()

    nlines = data.count('\n') - 2 # Last line is the footer
    del data
    fileobj.seek(start)
    data = pd.read_csv(fileobj, sep = "\t", nrows = nlines,
                       index_col = 0)
    new_index = pd.MultiIndex.from_tuples([(meta_data['sample'], i)\
                                          for i in data.columns], 
                                          names = ['sample', 'info'])
    data.columns = new_index
    return((data, meta_data))



def parseGPL(fileobj, size, chunksize = 100000):
    pass


def parseGSE(fname, chunksize = 100000):
    """
    """
    temp = gzip.open(fname, 'rb')
    sum_bytes = 0
    ent_all = []
    # A dictionary with keys being the platform and values being the gsms
    data_list = collections.defaultdict(list)
    data_coll = {}
    # Go through the file getting where the file starts
    while True:
        content = temp.read(chunksize)
        if content == "":
            break

        entities = re.compile(u'\\^(SAMPLE|PLATFORM)').finditer(content)
        matches = [(i.start() + sum_bytes, i.group()) for i in entities]
        ent_all.extend(matches)
        sum_bytes += chunksize

    #Read till end for the last one
    chr_to_read = [ent_all[i+1][0]-ent_all[i][0] for i\
                         in range(len(ent_all)-1)]
    chr_to_read.append(None)
    print("%i entities found" % len(ent_all))
    temp.rewind()

    for i in range(len(ent_all)):
        if ent_all[i][1] == '^PLATFORM':
            temp.seek(ent_all[i][0])
            parseGPL(temp, chr_to_read[i])
        elif ent_all[i][1] == '^SAMPLE':
            temp.seek(ent_all[i][0])
            data, meta = parseGSM(temp, chr_to_read[i])
            #data_list[meta['platform']].append(data)
            data_coll[meta['platform']] = data_coll.get(meta['platform'],
                                                        data).append(data)
            print(data_coll[meta['platform']].columns)



class GSM(object):
    """
    """
    def __init__(self, datatable, meta):
        self.data = datatable
        self.meta = meta

    def info(self):
        pass


class eSet(object):
    """ A python version of bioconductor's eset using panda's dataframes
    """
    def __init__(self):
        pass

    def exprs(self):
        pass
