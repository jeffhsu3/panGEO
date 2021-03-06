import urllib2, os, gzip, re, tempfile, cStringIO, StringIO, collections

import pandas as pd

import eSet


def _construct_geo_ftp(geo_id, amount = 'full'):
    """
    Construct an ftp address from the GEO ID

    """

    url_header = 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/'
    geo_type = geo_id[0:3]
    if geo_type == 'GSE':
        url_header += 'SOFT/by_series/%s/%s_family.soft.gz' % (geo_id,
                geo_id)
    elif geo_type == 'GPL':
        #:TODO finish this
        pass
    elif geo_type == 'GSM':
        #:TODO finish this
        pass
    elif geo_type == 'GDS':
        url_header += 'SOFT/GDS/%s.soft.gz' % geo_id
    print(url_header)
    return(url_header)


def _parse_meta(metalines):
    meta_info = collections.defaultdict(list)

    metalines = metalines.split("\n")
    m_pre = re.compile(u'!\w*\s=')
    metalines = [tuple(re.sub(u'!\w*?_','', i).split(" = "))\
                 for i in metalines if m_pre.search(i)]
    # Need to do this because some entries have multiple values
    for i in metalines:
        meta_info[i[0]].append(i[1])
    return(meta_info)


def getGEO(geo_id, out_file='/tmp/'):
    url_stream = urllib2.urlopen(_construct_geo_ftp(geo_id))
    # :TODO seems terribly inefficient to save and reload the gzip. Use
    # streams. Refactor to use a stream.
    out = os.path.join(out_file, geo_id)
    open(out, 'wb').write(url_stream.read())
    url_stream.close()
    return(parseGSE(out))


def parseGSM(fileobj, size, chunksize = 100000):
    """ Parse GSM.  Right now parse GSM and parse GSL are exactly the same.
    """
    metalines = ''
    while True:
        line = fileobj.readline()
        metalines += line
        if re.search('_table_begin', line):
            start = fileobj.tell()
            break
        else: pass
    seen = len(metalines)
    meta_data = _parse_meta(metalines)
    # Get the line counts to read
    if size:
        data = fileobj.read(size-seen)
    else:
        data = fileobj.read()

    nlines = data.count('\n') - 2 # Last line is the footer
    # :TODO since data is already read in, probably faster to parse this and
    # load it directly into a dataframe?
    del data
    fileobj.seek(start)
    # Can we make this faster?
    data = pd.read_csv(fileobj, sep = "\t", nrows = nlines,
                       index_col = 0)
    new_index = pd.MultiIndex.from_tuples([(meta_data['title'][0], i)\
                                          for i in data.columns],
                                          names = ['sample', 'info'])
    data.columns = new_index
    return((data, meta_data))


def parseGPL(fileobj, size, chunksize = 100000):
    """ Parses a GPL file.
    """
    metalines = ''
    while True:
        line = fileobj.readline()
        metalines += line
        if re.search('_table_begin', line):
            start = fileobj.tell()
            break
        else: pass
    seen = len(metalines)
    meta_data = _parse_meta(metalines)

    # Get the line counts to read
    if size:
        data = fileobj.read(size-seen)
    else:
        data = fileobj.read()

    meta_data = _parse_meta(metalines)
    nlines = data.count('\n') - 2 # Last line is the footer
    del data # Maybe easy to just parse and load this?
    fileobj.seek(start)
    gpl = pd.read_csv(fileobj, sep = "\t", nrows = nlines,
                      index_col = 0)
    return((gpl, meta_data))


def parseGSE_2(ibuffer):
    raise NotImplemented


def parseGSE(fname, chunksize = 100000):
    """
    Parses a GSE file retuns an eSet.

    GSE will often times have data from multiple platforms, so the data is
    returned as a dictionary with keys being the GPL id and the values being
    dataframes of all the data in that series from that one platform.

        >>> parseGSE()
    """

    temp = gzip.open(fname, 'rb')
    sum_bytes = 0
    ent_all = []
    # A dictionary with keys being the platform and values being the gsms
    data_list = collections.defaultdict(list)
    data_coll = {}
    annot_coll = {}
    meta_coll = {}
    # Go through the file getting where each of the samples starts
    while True:
        content = temp.read(chunksize)
        if content == "":
            break

        entities = re.compile(u'\\^(SAMPLE|PLATFORM)').finditer(content)
        matches = [(i.start() + sum_bytes, i.group()) for i in entities]
        ent_all.extend(matches)
        sum_bytes += chunksize

    chr_to_read = [ent_all[i+1][0]-ent_all[i][0] for i\
                         in range(len(ent_all)-1)]
    chr_to_read.append(None)
    temp.rewind()

    for i in range(len(ent_all)):
        if ent_all[i][1] == '^PLATFORM':
            temp.seek(ent_all[i][0])
            annot, gpl_meta = parseGPL(temp, chr_to_read[i])
            try:
                annot_coll = annot
            except KeyError:
                print('Key Error')
                pass
        elif ent_all[i][1] == '^SAMPLE':
            temp.seek(ent_all[i][0])
            data, meta = parseGSM(temp, chr_to_read[i])
            # :TODO meta needs to be turned into a dataframe?  Not sure since
            # there are duplicate rows for some.  However this is annoying have
            # to deal with them as a list.  Can't simply unlist entries that
            # only have 1 entry since not all fields are guranteed to have one
            # entry across all samples.
            # :TODO refactor this
            try:
                data_coll[meta['platform_id'][0]] =\
                        pd.concat([data_coll[meta['platform_id'][0]], data],
                                  axis=1)
            except KeyError:
                data_coll[meta['platform_id'][0]] = data
        else: pass

    temp.close()
    return(eSet.GSE(data_coll, annot_coll, meta_coll))
