import inputs_reader
from periodicity_tools import handle_periodicity

class InputData:
    """
        Class for input data: functions to read data, put them in xarrays, and manage indices of the domain (manage continuity around 0 and 180°) 
            - required arguments: 
            - optional keyword arguments: 
    """

    def __repr__(self):
        return str(self.__dict__.keys())

    def __init__(self, **kwargs):
        self.var = kwargs.pop('var', [])
        self.input_file = kwargs.pop('input_file', [])
        self.grid = kwargs.pop('grid', [])
        self.dico = kwargs.pop('dico', [])
        self.ds = []
        self.data = []

    def open_dataset(self):

        self.ds = xr.open_dataset(self.input_file)

    def select_area(self, crocogrd):

        londata = self.dico['lon' + self.grid]
        latdata = self.dico['lat' + self.grid]
        # manage 0-360 or -180-180
        if crocogrd.lonmin<0:
            # case 0-360 => to -180-180
            self.ds = d.assign_coords({londata: (((getattr(d, londata) + 180) % 360) - 180)}).sortby(londata)
            # or other solution: tmp = d.X.data; tmp[tmp>180] -= 360; d = d.assign_coords(X = tmp).sortby('X')
        else:
            # case -180-180 to 0-360
            self.ds = d.assign_coords({londata:  (getattr(d, londata) % 360)}).sortby(londata)

        # select sub-domain around the croco grid area
        self.ds = self.data.sel({latdata : slice(crocogrd.latmin-dl, crocogrd.latmax+dl),
                             londata : slice(crocogrd.lonmin-dl, crocogrd.lonmax+dl)})

    def get_data(self, crocogrd):
        self.open_dataset()
        self.select_area(crocogrd)
        self.data = getattr(self.ds, self.var)


