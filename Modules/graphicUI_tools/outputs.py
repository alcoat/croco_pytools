from traits.api import CArray, HasTraits

class Outputs(HasTraits):
    """
    Outputs object
    """
    lon_rho = CArray()
    lat_rho = CArray()
    lon_u = CArray()
    lat_u = CArray()
    lon_v = CArray()
    lat_v = CArray()
    h = CArray()
    hraw = CArray()
    pm = CArray()
    pn = CArray()
    angle = CArray()
    f = CArray()
    mask_rho = CArray()
