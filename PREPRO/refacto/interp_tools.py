import numpy as np
import scipy.interpolate as itp
import pyinterp.backends.xarray as pyxr
from progressbar import progressbar

def interp_horiz(myvar, crocogrd):
    '''
    Lon/Lat horizontal interpolation on croco grid 'r':
        - check the masked values within netcdf _FillValue, otherwise assume 0 is fill value
        - count the number of good/missing data: 
            - if no good data: put nan for this layer
            - if less than 10 good data: put the average value everywhere before interpolating
            - if 10 or more good data: compute ND nearest interpolator
        - finally perform the interpolation with bivariate pyinterp function  
    '''

    if "_FillValue" not in myvar.input_file.encoding:# If no FillValue in netcdf, assume 0 as value for the mask  
        myvar.data[myvar.data==0] = np.nan

    igood = np.where(np.isnan(myvar.data)==False)
    ibad  = np.where(np.isnan(myvar.data))
    NzGood = np.size(igood)
    Nbad = np.size(ibad)

    londata = myvar.dico['lon' + myvar.grid]
    latdata = myvar.dico['lat' + myvar.grid]

    if NzGood==0: tmpvar = np.nan
    else:
        if NzGood<10: # average the good values
            myvar.data = myvar.data.nanmean()
        else:
            spline = itp.NearestNDInterpolator( (getattr(myvar.ds, londata)[igood].ravel(),
                                                 getattr(myvar.ds, latdata)[igood].ravel() ),
                                                myvar.data[igood[0],igood[1]]
                                              )
            myvar.data[ibad] = spline(getattr(myvar.ds, londata)[ibad[0],ibad[1]],
                                      getattr(myvar.ds, latdata)[ibad[0],ibad[1]])

        val_interpolator = pyxr.Grid2D(myvar.data)
        tmpvar = val_interpolator.bivariate(coords=dict(lon=crocogrd.lon.flatten(),lat=crocogrd.lat.flatten()),num_threads=1).reshape(crocogrd.lon.shape)

    return tmpvar, NzGood

def interp_horiz_var3D(mydata, crocogrd):
    '''
    Interpolation on the croco grid or 3D variables: 
        - horizontal interpolation on each z level
        - concatenate all z levels 
        - add 2 levels (above and under) to avoid vertical extrap.
    '''

    vertvar = []
    Nz = np.shape(mydata.ds.depth)
    for k in progressbar(range(Nz), mydata.var + ': ', 40):#range(Nz)
        tmpvar, NzGood = interp_horiz(mydata.data[k,:,:], crocogrd)
        if Nzgood>Nzgoodmin: vertvar.append(tmpvar)
    vertvar = np.array(vertvar)

    Z = mydata.ds.depth[0:len(vertvar)]
    # put depth negative
    Z[Z>0] = Z[Z>0]*(-1)
    # flip the matrices upside down (Z[Nz]=surface)
    Z = np.flipud(np.concatenate(([100.], Z, [-10000.])))

    # Add a layer below the bottom and above the surface to avoid vertical extrapolations
    vertvar = np.concatenate((vertvar[0], vertvar, vertvar[-1]), axis=0)
    vertvar = np.flip(vertvar, axis=0)

    return vertvar, Z

