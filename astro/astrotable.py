from . import astrohelpers
from ..table import Table, __indent__
import numpy as np


class AstroTable(Table):
    """ Derived from the Table, this class add implementations of common astro tools especially conesearch """
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.ra_name, self.dec_name = self.__autoRADEC__()
        if (len(args) > 0):
            if isinstance(args[0], AstroTable):
                self.ra_name = args[0].ra_name
                self.dec_name = args[0].dec_name
        self.ra_name = kwargs.get('ra_name', self.ra_name)
        self.dec_name = kwargs.get('dec_name', self.dec_name)

    def __autoRADEC__(self):
        """ Tries to identify the columns containing RA and DEC coordinates """
        if 'ra' in self:
            ra_name = 'ra'
        elif 'RA' in self:
            ra_name = 'RA'
        else:
            ra_name = None
        if 'dec' in self:
            dec_name = 'dec'
        elif 'DEC' in self:
            dec_name = 'DEC'
        else:
            dec_name = None
        return ra_name, dec_name

    def setRA(self, val):
        """ Set the column that defines RA coordinates """
        assert(val in self), 'column name {} not found in the table'.format(val)
        self.ra_name = val

    def setDEC(self, val):
        """ Set the column that defines DEC coordinates """
        assert(val in self), 'column name {} not found in the table'.format(val)
        self.dec_name = val

    def getRA(self, degree=True):
        """ Returns RA, converted from hexa/sexa into degrees """
        if self.ra_name is None:
            return None
        if (not degree) or (self.dtype[self.ra_name].kind != 'S'):
            return self[self.ra_name]
        else:
            if (len(self[0][self.ra_name].split(':')) == 3):
                return np.asarray(astrohelpers.hms2deg(self[self.ra_name], delim=':'))
            elif (len(self[0][self.ra_name].split(' ')) == 3):
                return np.asarray(astrohelpers.hms2deg(self[self.ra_name], delim=' '))
            else:
                raise Exception('RA Format not understood')

    def getDEC(self, degree=True):
        """ Returns RA, converted from hexa/sexa into degrees """
        if self.dec_name is None:
            return None
        if (not degree) or (self.dtype[self.dec_name].kind != 'S'):
            return self[self.dec_name]
        else:
            if (len(self[0][self.dec_name].split(':')) == 3):
                return np.asarray(astrohelpers.dms2deg(self[self.dec_name], delim=':'))
            elif (len(self[0][self.dec_name].split(' ')) == 3):
                return np.asarray(astrohelpers.dms2deg(self[self.dec_name], delim=' '))
            else:
                raise Exception('RA Format not understood')

    def info(self):
        print self.header
        print "Table contains: %i row(s) in %i column(s)\n" % (self.nrows, self.ncols)
        if (self.ra_name is not None) & (self.dec_name is not None):
            print "Position coordinate columns: {}, {}\n".format(self.ra_name, self.dec_name)
        if self._aliases is not None:
                if len(self._aliases) > 0:
                        print "Table contains alias(es):"
                        for k, v in self._aliases.iteritems():
                                print '\t %s --> %s' % (k, v)
                        print ''
        fields = 'columns unit format description'.split()
        row    = [ (k, self.columns[k].unit, self.columns[k].format, self.columns[k].description) for k in self.keys() ]
        out    = __indent__([fields] + row, hasHeader=True, hasUnits=False, delim=' ')
        print out

    def coneSearch(self, ra, dec, r, outtype=0):
        """ Perform a cone search on a table
        INPUTS:
            ra0 	ndarray[ndim=1, dtype=float]	column name to use as RA source in degrees
            dec0	ndarray[ndim=1, dtype=float]	column name to use as DEC source in degrees
            ra		float                       	ra to look for (in degree)
            dec		float	                        ra to look for (in degree)
            r		float		                    distance in degrees
        KEYWORDS:
            outtype int                             0 -- minimal, indices of matching coordinates
                                                    1 -- indices and distances of matching coordinates
                                                    2 -- full, boolean filter and distances
        """

        assert( (self.ra_name is not None) & (self.dec_name is not None) ), 'Coordinate columns not set.'

        ra0  = self.getRA()
        dec0 = self.getDEC()
        return astrohelpers.conesearch(ra0, dec0, ra, dec, r, outtype=outtype)

    def selectWhere(self, fields, condition=None, condvars=None, cone=None, **kwargs):
        """ Read table data fulfilling the given `condition`.
            Only the rows fulfilling the `condition` are included in the result.
            conesearch is also possible trhough the keyword cone formatted as (ra, dec, r)
        """
        if (condition is None) & (cone is None):
            tab = self.__class__(self)

        if cone is not None:
            assert(len(cone) == 3), 'Expecting cone keywords as a triplet (ra, dec, r)'

        if condition is not None:
            tab = super(self.__class__, self).selectWhere(fields, condition, condvars, **kwargs)
            if cone is not None:
                ra, dec, r = cone
                b, d = tab.coneSearch(ra, dec, r, outtype=1)
                tab.data = tab.data[b]
                tab.add_column('separation', np.asarray(d), unit='degree')
                tab.header['COMMENT'] = 'SELECT %s FROM %s WHERE %s' % (fields, self.header['NAME'], 'distance from (%0.3f, %0.3f) <= %0.3f' % (ra, dec, r) )

        else:
            # only a cone search
            ra, dec, r = cone
            tab = self.__class__(self)
            ind, d = tab.coneSearch(ra, dec, r, outtype=1)

            if fields.count(',') > 0:
                _fields = fields.split(',')
            elif fields.count(' ') > 0:
                _fields = fields.split(' ')
            else:
                _fields = fields
            if _fields == '*':
                tab.data = tab.data[ind]
            else:
                tab.data = tab.data[tab.resolve_alias(_fields)][ind]
                names = tab.data.dtype.names
                #cleanup aliases and columns
                for k in self.keys():
                    if k not in names:
                        al = self.reverse_alias(k)
                        for alk in al:
                            self.delCol(alk)
                        tab.columns.pop(k)

            tab.add_column('separation', np.asarray(d), unit='degree')
            tab.header['COMMENT'] = 'SELECT %s FROM %s WHERE %s' % (fields, self.header['NAME'], 'distance from (%0.3f, %0.3f) <= %0.3f' % (ra, dec, r) )

        tab.setRA(self.ra_name)
        tab.setDEC(self.dec_name)
        return tab
