import pandas as pd
from scipy.stats import scoreatpercentile 

class Base_Array(pd.DataFrame):
    """ Generic class for holding array information from arrays all of the same
    tyep.  Subclass for specific arrays. 
    """
    
    def summary(self):
        means = self.apply(np.mean, 1)
        quart_one = self.apply(scoreatpercentile, 1)
        quart_three = self.apply(sco[good_probes, 0::2]eatpercentle, 1)
        summary = pd.DataFrame({'means': means, 'quartile_one': quart_one,
            'quartile_three': quart_three])
        self.summary = summary
        print(self.summary)

    def histogram(self, **kwargs):
        """ Plot a histogram of all the samples
        """
        hist_data = [np.histogram(self[x], **kwargs) for x in self]

        for i, j in hist_data:
            j = 0.5 * (j[1:] + j[:-1])
            plot(j, i, '-')


class Illumina_array(Array):
    """ Functions specific to the Illumina platform
    """

    def detection(self, threshold = 0.05, sample_number = None):
        """ Returns the dataframe with only rows that have significant
        quantification assay.  Used p-detection value.  
        """
        platforms = GSE.keys()

        if sample_number == None:
            sample_number = self.shape[1]/4

        good_probes = self.ix[:, 1::2].apply((lambda x: np.sum( x < threshold))
                \ sample_number)

        return self[good_probes, 0::2]

    def summary(self):
        means = self.ix[:,0::2].apply(np.mean, 1)


class eSet(object):
    """
    """
    def __init__(self, data, annot, meta = None, desc = None):
        self.data = data
        self.annot = annot
        self.meta = meta
        self.desc = desc
    '''
    def __setattr__(self, name, value):
        """ For interactive use
        """
        if name == '_data':
            super(eSet, self).__setattr__(name, value)
        else:
            try:
                existing = getattr(self, name)
    '''

    def __getattr__(self, name):
        if name in self.data.keys():
            return self.data[name]
        raise AttributeError("'%s' object has no attribute '%s'" %
                             (type(self).__name__, name))

    def exprs(self):
        pass


class GEO_base(object):
    """ Base class for all GEO objects including, GSE, GPL and GSM objects.
    """

    def __init__(self, data, meta):
        self.data = data
        self.meta = meta


class GSE(object):
    """ A class holding GSE objects
    """


    def __init__(self, gsms, gpls, meta = None):
        self.gsms = gsms
        self.gpls = gpls
        self.meta = meta

    '''

    def __setattr__(self, name, value):
        """ For interactive use
        """
        if name == '_data':
            super(eSet, self).__setattr__(name, value)
        else:
            try:
                existing = getattr(self, name)
                if name in self.gsms.keys():
                    self[name] = value
                else:
                    object.__setattr__(self, name, value)
            except (AttributeError, TypeError):
                object.__setattr__(self, name, value)


    def __getattr__(self, name):
        if name in self.gsms.keys():
            return self.gsms[name]
        raise AttributeError("'%s' object has no attribute '%s'" %
                             (type(self).__name__, name))
    '''
    def _get_samples(self):
        platforms = self.gsms.keys()
        self.samples = None


    def subset(self, platform, matrix_only = False):
        """ Subsets a GSE on a platform and returns a GSM object. Or simply an
        unannotated matrix.
        """
        return(GSM(self.gsms[platform], self.meta[platform]))


class GPL(GEO_base):
    """
    """

    def test_object(self):
        pass


class GSM(GEO_base):
    """
    """
    def __init__(self, gsm, meta):
        # GSM should be a dataframe
        self.gsm = gsm
        self.meta = meta

    def test_object(self):
        pass

