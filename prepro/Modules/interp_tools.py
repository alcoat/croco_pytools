import os
import sys
import time

import Cgrid_transformation_tools as grd_tools
import numpy as np
import pyinterp
import pyinterp.backends.xarray as pyxr
import pyinterp.fill as pyfi
import scipy.interpolate as itp
import xarray as xr
from progressbar import progressbar
from scipy.interpolate import griddata
from scipy.spatial import Delaunay
from sigmagrid_tools import ztosigma

#


def make_xarray(data, lon2D, lat2D):
    data_array = xr.DataArray(
        data,
        dims=["lat", "lon"],
        coords=[lat2D[:, 0].astype(np.float64), lon2D[0, :].astype(np.float64)],
    )
    data_array.lon.attrs["units"] = "degrees_east"
    data_array.lat.attrs["units"] = "degrees_north"
    return data_array


###############
# Get interpolation weight
###############
def get_tri_coef(X, Y, newX, newY, verbose=0):
    """
    Inputs:
        origin lon and lat 2d arrays (X,Y)
        child lon and lat 2d arrays (newX,newY)
    Ouputs:
        elem - pointers to 2d gridded data (at lonp,latp locations) from
            which the interpolation is computed (3 for each child point)
        coef - linear interpolation coefficients
    Use:
        To subsequently interpolate data from Fp to Fc, the following
        will work:      Fc  = sum(coef.*Fp(elem),3);  This line  should come in place of all
        griddata calls. Since it avoids repeated triangulations and tsearches (that are done
        with every call to griddata) it should be much faster.
    """

    Xp = np.array([X.ravel(), Y.ravel()]).T
    Xc = np.array([newX.ravel(), newY.ravel()]).T

    # Compute Delaunay triangulation
    if verbose == 1:
        tstart = tm.time()
    tri = Delaunay(Xp)
    if verbose == 1:
        print("Delaunay Triangulation", tm.time() - tstart)

    # Compute enclosing simplex and barycentric coordinate (similar to tsearchn in MATLAB)
    npts = Xc.shape[0]
    p = np.zeros((npts, 3))

    points = tri.points[tri.vertices[tri.find_simplex(Xc)]]
    if verbose == 1:
        tstart = tm.time()
    for i in progressbar(range(npts), "  Get_tri_coef: ", 40):

        if verbose == 1:
            print(np.float(i) / npts)

        if tri.find_simplex(Xc[i]) == -1:  # Point outside triangulation
            p[i, :] = p[i, :] * np.nan

        else:

            if verbose == 1:
                tstart = tm.time()
            A = np.append(np.ones((3, 1)), points[i], axis=1)
            if verbose == 1:
                print("append A", tm.time() - tstart)

            if verbose == 1:
                tstart = tm.time()
            B = np.append(1.0, Xc[i])
            if verbose == 1:
                print("append B", tm.time() - tstart)

            if verbose == 1:
                tstart = tm.time()
            p[i, :] = np.linalg.lstsq(A.T, B.T)[0]
            if verbose == 1:
                print("solve", tm.time() - tstart)

    if verbose == 1:
        print("Coef. computation 1", tm.time() - tstart)

    if verbose == 1:
        tstart = tm.time()
    elem = np.reshape(
        tri.vertices[tri.find_simplex(Xc)], (newX.shape[0], newY.shape[1], 3)
    )
    coef = np.reshape(p, (newX.shape[0], newY.shape[1], 3))
    if verbose == 1:
        print("Coef. computation 2", tm.time() - tstart)

    return (elem, coef)


#######################################


def get_delaunay_bry(lon_bry, lat_bry, inputfile, bdy):
    """
    This function computes the delaunay matrices for the interpolations
    at the boundaies

    Input:
      lon_bry      Longitudes of the boundary (vector).
      lat_bry      Latitudes of the boundary (vector).
      inputfile    netcdf structure poiting to the input file
      bdy          which boundary is done

    Output:
      LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,
      LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,
      LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry
    """
    comp_delaunay = 1

    #
    # get grid positions
    #
    LonU_bry, LatU_bry = getattr(inputfile, "lonU" + bdy), getattr(
        inputfile, "latU" + bdy
    )
    LonV_bry, LatV_bry = getattr(inputfile, "lonV" + bdy), getattr(
        inputfile, "latV" + bdy
    )
    LonT_bry, LatT_bry = getattr(inputfile, "lonT" + bdy), getattr(
        inputfile, "latT" + bdy
    )

    #
    # Get the 2D interpolation coefficients
    #

    if comp_delaunay == 1:

        print("\nCompute Delaunay triangulation from T-points to CROCO rho_points...")
        print("-------------------------------------------------------------------")
        [elemT_bry, coefT_bry] = get_tri_coef(LonT_bry, LatT_bry, lon_bry, lat_bry)
        coefnorm = np.sum(coefT_bry, axis=2)
        coefT_bry = coefT_bry / coefnorm[:, :, np.newaxis]

        print("\nCompute Delaunay triangulation from U-points to CROCO rho_points...")
        print("-------------------------------------------------------------------")
        [elemU_bry, coefU_bry] = get_tri_coef(LonU_bry, LatU_bry, lon_bry, lat_bry)
        coefnorm = np.sum(coefU_bry, axis=2)
        coefU_bry = coefU_bry / coefnorm[:, :, np.newaxis]

        print("\nCompute Delaunay triangulation from V-points to CROCO rho_points...")
        print("-------------------------------------------------------------------")
        [elemV_bry, coefV_bry] = get_tri_coef(LonV_bry, LatV_bry, lon_bry, lat_bry)
        coefnorm = np.sum(coefV_bry, axis=2)
        coefV_bry = coefV_bry / coefnorm[:, :, np.newaxis]

        # Save the Delaunay triangulation matrices
        np.savez(
            "coeffs_bry" + bdy + ".npz",
            coefT=coefT_bry,
            elemT=elemT_bry,
            coefU=coefU_bry,
            elemU=elemU_bry,
            coefV=coefV_bry,
            elemV=elemV_bry,
        )

    else:
        #
        # Load the Delaunay triangulation matrices
        #
        print("Load Delaunay triangulation...")
        data = np.load("coeffs_bry" + bdy + ".npz")
        coefT_bry = data["coefT_bry"]
        elemT_bry = data["elemT_bry"]
        coefU_bry = data["coefU_bry"]
        elemU_bry = data["elemU_bry"]
        coefV_bry = data["coefV_bry"]
        elemV_bry = data["elemV_bry"]

    print("Delaunay triangulation done")

    return (elemT_bry, coefT_bry, elemU_bry, coefU_bry, elemV_bry, coefV_bry)


####################################################################
def add2layers(vin):
    """
    Add a layer below the bottom and above the surface to avoid
    vertical extrapolations when doing a vertical interpolation

    Input:
      Vin    3 or 4D Variable

    Output:
      vout   Vin with 2 new layers above and below
    """
    if len(np.shape(vin)) == 3:
        [Nz, M, L] = np.shape(vin)
        vout = np.zeros((Nz + 2, M, L))

        vout[1:-1, :, :] = vin
        vout[0, :, :] = vin[0, :]
        vout[-1, :, :] = vin[-1, :]

    elif len(np.shape(vin)) == 4:
        [T, Nz, M, L] = np.shape(vin)
        vout = np.zeros((T, Nz + 2, M, L))

        vout[:, 1:-1, :, :] = vin
        vout[:, 0, :, :] = vin[:, 0, :]
        vout[:, -1, :, :] = vin[:, -1, :]
    return vout


######################################################################
def interp_tracers(inputfile, vname, k, crocogrd, dtmin, dtmax, prev=0, nxt=0, bdy=""):
    """
    Remove the missing values from a gridded 2D field
    and do an horizontal interpolation using pyinterp tools

    Inputs:
      inputfile     Input class containing useful tools (see input_class.py)
      vname         Variable's name to interpolate
      k             Depth index ( -1 when no depth component)
      crocogrd      CROCO grid class
      dtmin         Starting index of the time series
      dtmax         Ending index of the time series
      prev          If 1 duplicates the first index of the time series
      nxt           If 1 duplicates the last index of the time series
      bdy           Specify which boundary you are on. Leave empty if not

    Output:
      vout          Interpolation of the field time series on 2D CROCO grid
      Nzgood        Number of good points on the input level k
    """

    # 0: Read input data informations
    nc = inputfile.ncglo
    varinp = inputfile.var
    if vname in ["u", "ubar"]:
        Lon, Lat = getattr(inputfile, "lonU" + bdy), getattr(inputfile, "latU" + bdy)
    elif vname in ["v", "vbar"]:
        Lon, Lat = getattr(inputfile, "lonV" + bdy), getattr(inputfile, "latV" + bdy)
    else:
        Lon, Lat = getattr(inputfile, "lonT" + bdy), getattr(inputfile, "latT" + bdy)

    # 0: Precompute valid indices (igood) and bad indices (ibad)

    Vin_one = inputfile.var_periodicity(vname, 0, k, bdy=bdy)
    if (
        "_FillValue" not in inputfile.ncglo[vname].encoding
    ):  # If no FillValue in netcdf, assume 0 as value for the mask
        Vin_one[Vin_one == 0] = np.nan

    valid_mask = ~np.isnan(
        Vin_one[:, :]
    )  # Assume the mask is the same for all time steps
    igood = np.where(valid_mask)
    ibad = np.where(~valid_mask)
    NzGood = np.size(igood)
    Nbad = np.size(ibad)

    # 1: Read data
    if dtmin != dtmax:
        l = np.arange(dtmin, dtmax + 1)
    else:
        l = dtmin

    Vin = inputfile.var_periodicity(vname, l, k, bdy=bdy)
    if dtmin == dtmax:
        Vin = Vin[np.newaxis, :]

    if prev == 1:  # create overlap before when no data avail
        Vtmp = np.zeros([Vin.shape[0] + 1, Vin.shape[1], Vin.shape[2]])
        Vtmp[1:, :] = np.copy(Vin)
        Vtmp[0, :] = np.copy(Vin[0, :])
        Vin = Vtmp
        del Vtmp
    if nxt == 1:  # create overlap after when no data avail
        Vtmp = np.zeros([Vin.shape[0] + 1, Vin.shape[1], Vin.shape[2]])
        Vtmp[:-1, :] = np.copy(Vin)
        Vtmp[-1, :] = np.copy(Vin[-1, :])
        Vin = Vtmp
        del Vtmp

    # 2: Remove bad values (using nearest values)

    if bdy != "":
        if bdy == "S":
            bound = "south"
        elif bdy == "W":
            bound = "west"
        elif bdy == "E":
            bound = "east"
        elif bdy == "N":
            bound = "north"
        crocolon = getattr(crocogrd, "lon_" + bound)
        crocolat = getattr(crocogrd, "lat_" + bound)
    else:
        crocolon = crocogrd.lon
        crocolat = crocogrd.lat

    Vout = np.zeros([Vin.shape[0], crocolon.shape[0], crocolon.shape[1]])

    for tt in range(Vin.shape[0]):
        if NzGood == 0:
            #            print('\nWarning: no good data')
            Vin[:] = np.nan
        elif NzGood < 10:
            #            print('\nWarning: less than 10 good values')
            Vinfilled = np.copy(Vin[tt, :])
            Vinfilled[:] = np.nanmean(Vinfilled[igood])
            val_interpolator = pyinterp.backends.xarray.Grid2D(
                make_xarray(Vinfilled, Lon, Lat)
            )
        else:
            spline = itp.NearestNDInterpolator(
                (Lon[igood].ravel(), Lat[igood].ravel()), Vin[tt, igood[0], igood[1]]
            )
            Vinfilled = np.copy(np.squeeze(Vin[tt, :]))
            Vinfilled[ibad] = spline(Lon[ibad[0], ibad[1]], Lat[ibad[0], ibad[1]])
            val_interpolator = pyinterp.backends.xarray.Grid2D(
                make_xarray(Vinfilled, Lon, Lat)
            )
        # 3: 2D interpolation
        if NzGood == 0:
            Vout[tt, :] = np.nan
        else:
            Vout[tt, :] = val_interpolator.bicubic(
                coords=dict(lon=crocolon.flatten(), lat=crocolat.flatten()),
                num_threads=1,
            ).reshape(crocolon.shape)
    return Vout, NzGood


#######################################################


def interp(
    inputfile, vname, Nzgoodmin, z_rho, crocogrd, dtmin, dtmax, prev=0, nxt=0, bdy=""
):
    """
    Do a full interpolation of a 3-4d variable from z-grid data to a CROCO sigma grid
    1 - Horizontal interpolation on each z-levels
    2 - Vertical Interpolation from z to CROCO sigma levels

    Inputs:
      inputfile     Input class containing useful tools (see input_class.py)
      vname         Variable's name to interpolate
      Nzgoodmin     Number of value to consider a z-level fine to be used
      z_rho         CROCO vertical levels
      crocogrd      CROCO grid class
      dtmin         Starting index of the time series
      dtmax         Ending index of the time series
      prev          If 1 duplicates the first index of the time series
      nxt           If 1 duplicates the last index of the time series
      bdy           Specify which boundary you are on. Leave empty if not

    Output:
      vout          3-4D interpolation of vname variable
    """

    if np.ndim(z_rho) == 3:
        z_rho = z_rho[np.newaxis, :]
    [T, N, M, L] = np.shape(z_rho)
    depth = inputfile.depth
    [Nz] = np.shape(depth)

    print("Horizontal interpolation over z levels")
    t4d = np.zeros((T, Nz, M, L))
    kgood = -1
    for k in progressbar(range(Nz), vname + ": ", 40):  # range(Nz)
        (t3d, Nzgood) = interp_tracers(
            inputfile, vname, k, crocogrd, dtmin, dtmax, prev, nxt, bdy=bdy
        )
        if Nzgood > Nzgoodmin:
            kgood = kgood + 1
            t4d[:, kgood, :, :] = t3d

    t4d = t4d[:, 0:kgood, :, :]
    depth = depth[0:kgood]

    [Nz] = np.shape(depth)
    if depth[0] < 0:
        Z = depth
    else:
        Z = -depth

    #  Vertical interpolation
    print("Vertical interpolation")

    # Add a layer below the bottom and above the surface to avoid vertical extrapolations
    # and flip the matrices upside down (Z[Nz]=surface)

    Z = np.flipud(np.concatenate(([100.0], Z, [-10000.0])))
    [Nz] = np.shape(Z)
    t4d = np.flip(add2layers(t4d), axis=1)

    # Do the vertical interpolations
    vout = ztosigma(t4d, Z, z_rho)
    return vout


####################################################################


def interp_uv(
    inputfile,
    Nzgoodmin,
    z_rho,
    cosa,
    sina,
    crocogrd,
    dtmin,
    dtmax,
    prev=0,
    nxt=0,
    bdy="",
):
    """
    Same as interp but for horizontal velocities u,v

    Inputs:
      inputfile     Input class containing useful tools (see input_class.py)
      Nzgoodmin     Number of value to consider a z-level fine to be used
      z_rho         CROCO vertical levels
      cosa          Cosine value of grid angle
      sina          Sinus value of grid angle
      crocogrd      CROCO grid class
      dtmin         Starting index of the time series
      dtmax         Ending index of the time series
      prev          If 1 duplicates the first index of the time series
      nxt           If 1 duplicates the last index of the time series
      bdy           Specify which boundary you are on. Leave empty if not

    Outputs:
       uout         U velocity on 4D (Time,CROCO-Ugrid)
       vout         V velocity on 4D (Time,CROCO-Vgrid)
       ubar         Integrated U velocity on (Time,CROCO-Ugrid)
       vbar         Integrated V velocity on (Time,CROCO-Vgrid)
    """

    if np.ndim(z_rho) == 3:
        z_rho = z_rho[np.newaxis, :]
    [T, N, M, L] = np.shape(z_rho)

    cosa3d = np.tile(cosa, (T, 1, 1))
    sina3d = np.tile(sina, (T, 1, 1))

    depth = inputfile.depth
    dz = np.gradient(depth)
    [Nz] = np.shape(depth)

    print("Horizontal interpolation of u and v over z levels")
    u4d = np.zeros((T, Nz, M, L - 1))
    v4d = np.zeros((T, Nz, M - 1, L))
    kgood = -1
    for k in progressbar(range(Nz), " uv : ", 40):
        (u3d, Nzgood_u) = interp_tracers(
            inputfile, "u", k, crocogrd, dtmin, dtmax, prev, nxt, bdy=bdy
        )
        (v3d, Nzgood_v) = interp_tracers(
            inputfile, "v", k, crocogrd, dtmin, dtmax, prev, nxt, bdy=bdy
        )
        Nzgood = np.min((Nzgood_u, Nzgood_v))
        if Nzgood > Nzgoodmin:
            kgood = kgood + 1
            #
            # Rotation and put to u-points and v-points
            #
            u4d[:, kgood, :, :] = grd_tools.rho2u(u3d * cosa3d + v3d * sina3d)
            v4d[:, kgood, :, :] = grd_tools.rho2v(v3d * cosa3d - u3d * sina3d)

    u4d = u4d[:, 0:kgood, :, :]
    v4d = v4d[:, 0:kgood, :, :]
    depth = depth[0:kgood]

    [Nz] = np.shape(depth)
    if depth[0] < 0:
        Z = depth
    else:
        Z = -depth

    # --------------------------------------------------
    #  Vertical interpolation
    # ----------------------------------------------------
    #
    print("Vertical interpolation")
    #
    # Add a layer below the bottom and above the surface to avoid vertical extrapolations
    # and flip the matrices upside down (Z[Nz]=surface)
    #
    Z = np.flipud(np.concatenate(([100.0], Z, [-10000.0])))
    [Nz] = np.shape(Z)
    u4d = np.flip(add2layers(u4d), axis=1)
    v4d = np.flip(add2layers(v4d), axis=1)
    #
    # Do the vertical interpolations
    #
    uout = ztosigma(u4d, Z, grd_tools.rho2u(z_rho))
    vout = ztosigma(v4d, Z, grd_tools.rho2v(z_rho))

    return uout, vout


####################################################################


def interp_tides(
    inputfile, vname, k, crocogrd, dtmin, dtmax, input_type, prev=0, nxt=0, bdy=""
):
    """
    Compute Extrapolate data with pyinterp gauss_seidel and perform horizontal interpolation for tides

    Inputs:
      inputfile     Input class containing useful tools (see input_class.py)
      vname         Variable's name to interpolate
      Nzgoodmin     Number of value to consider a z-level fine to be used
      z_rho         CROCO vertical levels
      coef          Coef from Delaunay triangulation
      elem          Elem from Delaunay triangulation
      dtmin         Starting index of the time series
      dtmax         Ending index of the time series
      prev          If 1 duplicates the first index of the time series
      nxt           If 1 duplicates the last index of the time series
      bdy           Specify which boundary you are on. Leave empty if not

    Output:
      vout          4D interpolation of vname variable (Time+ CROCO-grid)
    """
    # 0: Read input data informations
    nc = inputfile.ncglo
    varinp = inputfile.var

    if vname == "H":
        var = ["ssh"]
    elif vname == "cur":
        var = ["u", "v"]

    for field in var:
        if field == "u":
            Lon, Lat = getattr(inputfile, "lonU" + bdy), getattr(
                inputfile, "latU" + bdy
            )
        elif field == "v":
            Lon, Lat = getattr(inputfile, "lonV" + bdy), getattr(
                inputfile, "latV" + bdy
            )
        else:
            Lon, Lat = getattr(inputfile, "lonT" + bdy), getattr(
                inputfile, "latT" + bdy
            )

        # 1: Read data
        if dtmin != dtmax:
            l = np.arange(dtmin, dtmax + 1)
        else:
            l = dtmin

        part1 = inputfile.var_periodicity(field + "_part1", l, k, bdy=bdy)
        part2 = inputfile.var_periodicity(field + "_part2", l, k, bdy=bdy)

        if input_type == "Amp_Phase":
            rpart, impart = part1 * np.cos(part2), part1 * np.sin(part2)
        else:
            rpart, impart = part1, part2

        if dtmin == dtmax:
            rpart = rpart[np.newaxis, :]
            impart = impart[np.newaxis, :]

        # 2: Remove bad values (using nearest values)
        #        if "_FillValue" not in inputfile.ncglo[field+'_part1'].encoding:# If no FillValue in netcdf, assume 0 as value for the mask
        #            rpart[rpart==0]=np.nan
        #            impart[impart==0]=np.nan

        if bdy != "":
            if bdy == "S":
                bound = "south"
            elif bdy == "W":
                bound = "west"
            elif bdy == "E":
                bound = "east"
            elif bdy == "N":
                bound = "north"
            crocolon = getattr(crocogrd, "lon_" + bound)
            crocolat = getattr(crocogrd, "lat_" + bound)
        else:
            crocolon = crocogrd.lon
            crocolat = crocogrd.lat

        Rout = np.zeros([rpart.shape[0], crocolon.shape[0], crocolon.shape[1]])
        Iout = np.zeros([impart.shape[0], crocolon.shape[0], crocolon.shape[1]])

        for tt in range(rpart.shape[0]):
            # assume rpart and impart have the same mask
            igood = np.where(np.isnan(rpart[tt, :]) == False)
            ibad = np.where(np.isnan(rpart[tt, :]))

            NzGood = np.size(igood)
            Nbad = np.size(ibad)

            if NzGood == 0:
                rpart[:] = np.nan
                impart[:] = np.nan
            elif NzGood < 10:
                rpartfilled = np.copy(rpart[tt, :])
                rpartfilled[:] = np.nanmean(rpartfilled[igood])
                r_val_interpolator = pyinterp.backends.xarray.Grid2D(
                    make_xarray(rpartfilled, Lon, Lat)
                )
                impartfilled = np.copy(impart[tt, :])
                impartfilled[:] = np.nanmean(impartfilled[igood])
                i_val_interpolator = pyinterp.backends.xarray.Grid2D(
                    make_xarray(impartfilled, Lon, Lat)
                )
            else:
                r_spline = itp.NearestNDInterpolator(
                    (Lon[igood].ravel(), Lat[igood].ravel()),
                    rpart[tt, igood[0], igood[1]],
                )
                rpartfilled = np.copy(np.squeeze(rpart[tt, :]))
                rpartfilled[ibad] = r_spline(
                    Lon[ibad[0], ibad[1]], Lat[ibad[0], ibad[1]]
                )
                r_val_interpolator = pyinterp.backends.xarray.Grid2D(
                    make_xarray(rpartfilled, Lon, Lat)
                )
                i_spline = itp.NearestNDInterpolator(
                    (Lon[igood].ravel(), Lat[igood].ravel()),
                    impart[tt, igood[0], igood[1]],
                )
                impartfilled = np.copy(np.squeeze(impart[tt, :]))
                impartfilled[ibad] = i_spline(
                    Lon[ibad[0], ibad[1]], Lat[ibad[0], ibad[1]]
                )
                i_val_interpolator = pyinterp.backends.xarray.Grid2D(
                    make_xarray(impartfilled, Lon, Lat)
                )
            # 3: 2D interpolation
            if NzGood == 0:
                Rout[tt, :] = np.nan
                Iout[tt, :] = np.nan
            else:
                Rout[tt, :] = r_val_interpolator.bicubic(
                    coords=dict(lon=crocolon.flatten(), lat=crocolat.flatten()),
                    num_threads=1,
                ).reshape(crocolon.shape)
                Iout[tt, :] = i_val_interpolator.bicubic(
                    coords=dict(lon=crocolon.flatten(), lat=crocolat.flatten()),
                    num_threads=1,
                ).reshape(crocolon.shape)

        Vout = Rout + 1j * Iout

        if field == "u":
            Uout = np.copy(Vout)

    if vname == "H":
        return Vout, NzGood

    elif vname == "cur":
        return Uout, Vout, NzGood


####################################################################


def compute_uvbar_ogcm(
    inputfile, cosa, sina, crocogrd, dtmin, dtmax, prev=0, nxt=0, bdy=""
):
    """
    Compute ubar and vbar from OGCM
    Remove the missing values (extrapolation)
    Do the horizontal interpolation using pyinterp tools to the croco grid
    Do the rotation and the shift on u- and v- croco grid

    Inputs:
    inputfile     Input class containing useful tools (see input_class.py)
    k             Depth index ( -1 when no depth component)
    crocogrd      CROCO grid class
    dtmin         Starting index of the time series
    dtmax         Ending index of the time series
    prev          If 1 duplicates the first index of the time series
    nxt           If 1 duplicates the last index of the time series
    bdy           Specify which boundary you are on. Leave empty if not

    Output:
    ubar_ogcm, vbar_ogcm, ubar_croco_rotated_u, vbar_croco_rotated_v

            ubar_ogcm, vbar_ogcm :
                        ubar and vbar on the OGCM grid computed
                        from the reanalysis (lon/lat grid)

                    ubar_croco_rotated_u, vbar_croco_rotated_v :
                        ubar  and vbar computed from the OGCM
                        and interpolated on the croco grid,
                        rotated, at U- and V- grid

    """
    # 0: Read input data informations
    nc = inputfile.ncglo
    varinp = inputfile.var

    Lon_ogcm = getattr(inputfile, "lonT" + bdy)
    Lat_ogcm = getattr(inputfile, "latT" + bdy)

    # 1: Read data
    if dtmin != dtmax:
        l = np.arange(dtmin, dtmax + 1)
    else:
        l = dtmin

    zeta_ogcm = inputfile.var_periodicity("ssh", l, -1, bdy=bdy)

    if np.ndim(zeta_ogcm) == 2:  # add a dimension of 1 in case of ini
        zeta_ogcm = zeta_ogcm[np.newaxis, :]

    [T, M, L] = np.shape(zeta_ogcm)

    depth = inputfile.depth
    [Nz] = np.shape(depth)

    print("Compute ubar and vbar from OGCM over z-levels and interpolate to croco grid")
    u_ogcm = np.zeros((T, Nz, M, L))
    v_ogcm = np.zeros((T, Nz, M, L))
    dz = np.gradient(depth)

    # for k in range(Nz):
    for k in progressbar(range(Nz), " - reading/extracting u/v from OGCM : ", 40):
        # print('K=',k)
        u_ogcm_tmp = inputfile.var_periodicity("u", l, k, bdy=bdy)
        v_ogcm_tmp = inputfile.var_periodicity("v", l, k, bdy=bdy)

        u_ogcm[:, k, :, :] = u_ogcm_tmp
        v_ogcm[:, k, :, :] = v_ogcm_tmp

    # 1: Compute the barotropic velocities from reanalysis
    maskuv3d = np.where(np.isnan(u_ogcm), 0, 1)
    u_ogcm[np.isnan(u_ogcm)] = 0
    v_ogcm[np.isnan(v_ogcm)] = 0

    du = np.zeros((T, M, L))
    dv = np.zeros((T, M, L))
    zuv = du

    # for k in range(Nz):
    for k in progressbar(range(Nz), " - computing ubar/vbar from OGCM : ", 40):
        # print('K=',k)
        du = du + dz[k] * u_ogcm[:, k, :, :]
        dv = dv + dz[k] * v_ogcm[:, k, :, :]
        zuv = zuv + dz[k] * maskuv3d[:, k, :, :]

    du[zuv == 0] = np.nan
    dv[zuv == 0] = np.nan
    zuv[zuv == 0] = np.nan

    ubar_ogcm = du / zuv
    vbar_ogcm = dv / zuv

    # 0: Precompute valid indices (igood) and bad indices (ibad)
    ubar_ogcm_one = ubar_ogcm[0, :, :]
    vbar_ogcm_one = vbar_ogcm[0, :, :]
    igood = np.where(np.isnan(ubar_ogcm_one[:]) == False)
    ibad = np.where(np.isnan(ubar_ogcm_one[:]))
    NzGood = np.size(igood)
    Nbad = np.size(ibad)

    # 0: Precompute the interpolator spline only once
    splineU = itp.NearestNDInterpolator(
        (Lon_ogcm[igood].ravel(), Lat_ogcm[igood].ravel()),
        ubar_ogcm_one[igood[0], igood[1]],
    )

    splineV = itp.NearestNDInterpolator(
        (Lon_ogcm[igood].ravel(), Lat_ogcm[igood].ravel()),
        vbar_ogcm_one[igood[0], igood[1]],
    )

    if prev == 1:  # create overlap before when no data avail
        vartmp = np.zeros(
            [ubar_ogcm.shape[0] + 1, ubar_ogcm.shape[1], ubar_ogcm.shape[2]]
        )
        vartmp[1:, :] = np.copy(ubar_ogcm)
        vartmp[0, :] = np.copy(ubar_ogcm[0, :])
        ubar_ogcm = vartmp
        del vartmp

        vartmp = np.zeros(
            [vbar_ogcm.shape[0] + 1, vbar_ogcm.shape[1], vbar_ogcm.shape[2]]
        )
        vartmp[1:, :] = np.copy(vbar_ogcm)
        vartmp[0, :] = np.copy(vbar_ogcm[0, :])
        vbar_ogcm = vartmp
        del vartmp

    if nxt == 1:  # create overlap after when no data avail
        vartmp = np.zeros(
            [ubar_ogcm.shape[0] + 1, ubar_ogcm.shape[1], ubar_ogcm.shape[2]]
        )
        vartmp[:-1, :] = np.copy(ubar_ogcm)
        vartmp[-1, :] = np.copy(ubar_ogcm[-1, :])
        ubar_ogcm = vartmp
        del vartmp

        vartmp = np.zeros(
            [vbar_ogcm.shape[0] + 1, vbar_ogcm.shape[1], vbar_ogcm.shape[2]]
        )
        vartmp[:-1, :] = np.copy(vbar_ogcm)
        vartmp[-1, :] = np.copy(vbar_ogcm[-1, :])
        vbar_ogcm = vartmp
        del vartmp

    cosa3d = np.tile(cosa, (vbar_ogcm.shape[0], 1, 1))
    sina3d = np.tile(sina, (vbar_ogcm.shape[0], 1, 1))

    # print(' - horizontal interpolation to croco grid')
    # 2: Interp of ubar and vbar from reanalysis (ubar_ogcm and vbar_ogcm) to crocogrid
    if bdy != "":
        if bdy == "S":
            bound = "south"
        elif bdy == "W":
            bound = "west"
        elif bdy == "E":
            bound = "east"
        elif bdy == "N":
            bound = "north"

        crocolon = getattr(crocogrd, "lon_" + bound)
        crocolat = getattr(crocogrd, "lat_" + bound)
        crocolonu = getattr(crocogrd, "lonu_" + bound)
        crocolatu = getattr(crocogrd, "latu_" + bound)
        crocolonv = getattr(crocogrd, "lonv_" + bound)
        crocolatv = getattr(crocogrd, "latv_" + bound)
    else:
        crocolon = crocogrd.lon
        crocolat = crocogrd.lat
        crocolonu = crocogrd.lonu
        crocolatu = crocogrd.latu
        crocolonv = crocogrd.lonv
        crocolatv = crocogrd.latv

    ubar_croco_r = np.zeros([ubar_ogcm.shape[0], crocolon.shape[0], crocolon.shape[1]])
    vbar_croco_r = np.zeros([vbar_ogcm.shape[0], crocolon.shape[0], crocolon.shape[1]])
    ubar_croco_rotated_u = np.zeros(
        [ubar_ogcm.shape[0], crocolonu.shape[0], crocolonu.shape[1]]
    )
    vbar_croco_rotated_v = np.zeros(
        [vbar_ogcm.shape[0], crocolonv.shape[0], crocolonv.shape[1]]
    )

    # 2.1 2D extrapolation of ubar_ogcm and vbar_ogcm
    # for tt in range(ubar_croco_r.shape[0]):
    for tt in progressbar(
        range(ubar_croco_r.shape[0]), " -  horizontal interpolation to croco grid", 40
    ):
        #  ==> for ubar_ogcm
        ubar_ogcm_filled = np.copy(np.squeeze(ubar_ogcm[tt, :]))
        ubar_ogcm_filled[ibad] = splineU(
            Lon_ogcm[ibad[0], ibad[1]], Lat_ogcm[ibad[0], ibad[1]]
        )
        val_interpolatorU = pyinterp.backends.xarray.Grid2D(
            make_xarray(ubar_ogcm_filled, Lon_ogcm, Lat_ogcm)
        )

        #  ==>  for vbar_ogcm
        vbar_ogcm_filled = np.copy(np.squeeze(vbar_ogcm[tt, :]))
        vbar_ogcm_filled[ibad] = splineV(
            Lon_ogcm[ibad[0], ibad[1]], Lat_ogcm[ibad[0], ibad[1]]
        )
        val_interpolatorV = pyinterp.backends.xarray.Grid2D(
            make_xarray(vbar_ogcm_filled, Lon_ogcm, Lat_ogcm)
        )

        # 2.2 Interpolation to croco grid
        ubar_croco_r[tt, :] = val_interpolatorU.bicubic(
            coords=dict(lon=crocolon.flatten(), lat=crocolat.flatten()), num_threads=1
        ).reshape(crocolon.shape)
        vbar_croco_r[tt, :] = val_interpolatorV.bicubic(
            coords=dict(lon=crocolon.flatten(), lat=crocolat.flatten()), num_threads=1
        ).reshape(crocolon.shape)
        # 2.3 Rotation and shift on u- and v- grid
        ubar_croco_rotated_u[tt, :] = grd_tools.rho2u(
            ubar_croco_r[tt, :] * cosa3d[tt, :] + vbar_croco_r[tt, :] * sina3d[tt, :]
        )
        vbar_croco_rotated_v[tt, :] = grd_tools.rho2v(
            vbar_croco_r[tt, :] * cosa3d[tt, :] - ubar_croco_r[tt, :] * sina3d[tt, :]
        )

    return ubar_ogcm, vbar_ogcm, ubar_croco_rotated_u, vbar_croco_rotated_v
