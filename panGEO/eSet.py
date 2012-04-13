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

