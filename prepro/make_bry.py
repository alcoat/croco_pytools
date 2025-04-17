#!/usr/bin/env python

"""
===========================================================================
Further Information:
  http://www.croco-ocean.org

This file is part of CROCOTOOLS

Create a CROCO bounday files
In the current state the script can handle:
    - mercator (glorys)
    - soda
    - eccov4  (V4.3)

To add a new dataset you just have to go in Readers/inputs_readers.py and
create a dico with correspondance between croco_var(ssh,u,v,temp,salt)
and dataset variable names.

The periodicity of a dataset (if it is -180/180 or 0/360) is handle by the
script. Just download the dataset and it should be fine

The script works as follow:
    - reads croco_grd
    - reads input data and restricts the area of spatial coverage
    - creates croco_bry.nc
    - Loop on open boundaries with:
          check if there is data before or after for continuity, if not duplicate first or last
          loop on var with:
              * horizontal interpolation
              * vertical interpolation
    - Writes data in netcdf

===========================================================================
"""
__author__ = "Mathieu Le Corre"
__email__ = "mathieu.le.corre@shom.fr"
__date__ = "2022-09"
__license__ = "GPL3"

# --- Dependencies ---------------------------------------------------------

import glob as glob
import pathlib
import sys

import Cgrid_transformation_tools as grd_tools
import croco_class as Croco
import ibc_class as Inp
import interp_tools
import netCDF4 as netcdf
import numpy as np
import pandas
import pylab as plt
import sigmagrid_tools as sig_tools
import xarray as xr
import yaml
from dateutil.relativedelta import relativedelta

sys.path.append("./Modules/")
sys.path.append("./Readers/")


def run_make_bry():
    with open("make_bry_def.yml", "r", encoding="utf8") as infile:
        bry_def = yaml.safe_load(infile)

    origindate = pandas.Timestamp(bry_def["origindate"])
    # Put origin date to the right format
    day_zero_num = origindate.to_pydatetime()
    day_zero_num = plt.date2num(day_zero_num)

    # Put start and end date to the right format
    start_date = pandas.Timestamp(bry_def["begindate"])

    dtstrdt = start_date.to_pydatetime()

    end_date = pandas.Timestamp(bry_def["enddate"])
    dtenddt = end_date.to_pydatetime()

    dtstr, dtend = plt.date2num(dtstrdt), plt.date2num(dtenddt)

    # --- Load croco_grd --------------------------------------------------
    croco_dir = pathlib.Path(bry_def["croco_dir"])
    grid_pathname = croco_dir / bry_def["croco_grd"]
    crocogrd = Croco.CROCO_grd(grid_pathname, bry_def["sigma_params"])

    # --- Initialize boundary vars ----------------------------------------

    crocogrd.WEST_grid()
    crocogrd.EAST_grid()
    crocogrd.SOUTH_grid()
    crocogrd.NORTH_grid()

    # --- Initialize input data class -------------------------------------

    inpdat = Inp.getdata(
        bry_def["inputdata"],
        bry_def["input_file"],
        crocogrd,
        bry_def["multi_files"],
        bry_def["tracers"],
        bdy=[bry_def["obc_dict"], bry_def["cycle_bry"]],
    )

    # --- Work on date format for the loop in time ------------------------
    output_file_format = bry_def["output_file_format"].upper()
    startloc = start_date
    # startloc = pandas.offsets.MonthBegin().rollback(start_date)
    if output_file_format == "DAILY":
        endloc = startloc + pandas.Timedelta(1, unit="day")
    elif output_file_format == "MONTHLY":
        startloc = pandas.offsets.MonthBegin().rollback(start_date)
        endloc = startloc + pandas.tseries.offsets.DateOffset(months=1)
    elif output_file_format == "YEARLY":
        startloc = pandas.offsets.YearBegin().rollback(start_date)
        endloc = startloc + pandas.tseries.offsets.DateOffset(years=1)
        # if plt.datetime.datetime.strptime(
        #    str(startloc), "%Y-%m-%d %H:%M:%S"
        # ).year == int(end_date.year):
        #    endloc = plt.num2date(dtend).replace(tzinfo=None)
        # else:
        #    endloc = plt.datetime.datetime(int(start_date[:4]), 12, 31, 12)

    elif output_file_format == "FULL":
        endloc = end_date  # plt.num2date(dtend).replace(tzinfo=None)
    else:
        print(
            f'\n Output file format "{output_file_format}" is not setup. Pease change it to DAILY, MONTHLY, YEARLY or FULL'
        )
        sys.exit()

    # --- Start time loop loop in time ------------------------------------

    while endloc <= end_date:

        # Load full time dataset
        time = pandas.DatetimeIndex(inpdat.ncglo["time"].values)
        # find index for the time range
        ind = np.where((time > startloc) & (time <= endloc))

        if len(ind[0]) == 0:
            print("\nData is missing for range %s to %s" % (startloc, endloc))
            sys.exit()

        [dtmin, dtmax] = np.min(ind), np.max(ind)

        # create file
        tmp_date = startloc
        # file name depending on format chosen
        bdy_filename = croco_dir / bry_def["bry_filename"].format(
            date_b=startloc, date_e=endloc
        )

        Croco.CROCO.create_bry_nc(
            None,
            bdy_filename,
            crocogrd,
            bry_def["obc_dict"],
            bry_def["cycle_bry"],
            tracers=bry_def["tracers"],
        )
        #
        print("\n-----------------------------------")
        print(
            " Processing Year %s - Month %02i - Day : %02i"
            % (startloc.year, startloc.month, startloc.day)
        )
        print(
            "         to Year %s - Month %02i - Day : %02i"
            % (endloc.year, endloc.month, endloc.day)
        )

        # Check if there is at least 1 point by month when doing Yearly or full
        if bry_def["output_file_format"].upper() != "MONTHLY":
            tmp_str_date = startloc
            tmp_end_date = tmp_str_date + relativedelta(months=1, days=-1, hours=12)
            while plt.date2num(tmp_end_date) <= plt.date2num(endloc):
                ind_tmp = np.where((time > tmp_str_date) & (time <= tmp_end_date))
                tmp_date = plt.datetime.datetime.strptime(
                    str(tmp_str_date), "%Y-%m-%d %H:%M:%S"
                )
                if len(ind_tmp[0]) == 0:
                    print(
                        "\nLacking %s data for  Y%s - M%02i"
                        % (bry_def["inputdata"], tmp_date.year, tmp_date.month)
                    )
                    sys.exit()
                tmp_str_date = tmp_end_date + relativedelta(days=1, hours=-12)
                tmp_end_date = tmp_str_date + relativedelta(months=1, days=-1, hours=12)

            del tmp_str_date, tmp_end_date

        # --- Check if data availabla for the surrounded months -----------
        print("-----------------------------------")
        tmp_str_date = startloc
        tmp_end_date = endloc
        prev_month_str = tmp_str_date + relativedelta(months=-1)
        prev_month_end = tmp_str_date + relativedelta(days=-1, hours=12)
        ind_prev = np.where((time > prev_month_str) & (time < prev_month_end))
        if len(ind_prev[0]) == 0:
            print("   No data for the previous month: using current month")
            prev = 1
        else:
            prev = 0
            dtmin = (
                dtmin - 1
            )  # create overlap before (in this case it takes the previous month)
        del prev_month_str, prev_month_end, ind_prev
        #
        next_month_str = tmp_end_date + relativedelta(days=1)
        next_month_end = next_month_str + relativedelta(months=1, days=-1, hours=12)
        ind_next = np.where((time > next_month_str) & (time < next_month_end))
        if len(ind_next[0]) == 0:
            print("   No data for the next month: using current month")
            nxt = 1
        else:
            nxt = 0
            dtmax = (
                dtmax + 1
            )  # create overlap after (in this case it takes the next month)

        del next_month_str, next_month_end, ind_next, tmp_str_date, tmp_end_date

        if time[dtmin : dtmax + 1].std() >= pandas.Timedelta(
            5, unit="days"
        ):  # Abnormal distribution of days
            Question = (
                input(
                    "Abnormal distribution of days (variance to high) \
                    \nThis may be due to the use of different temproral resolution dataset.\
                    \n Do you want to proceed?: y,[n] "
                )
                or "no"
            )
            if Question.lower() == ("n") or Question.lower() == ("no"):
                print("Aborting")
                sys.exit()

        ## --- Handle bry_time --------------------------------------------

        bry_time = (time[dtmin : dtmax + 1] - origindate).total_seconds() / 86400.0
        if prev == 1 and len(bry_time) == 1:
            prev_time = plt.date2num(
                plt.num2date(bry_time[0]) + relativedelta(days=-30)
            )
            bry_time = np.append(prev_time, bry_time)
            del prev_time
        elif prev == 1 and len(bry_time) > 1:
            date_dt = np.gradient(bry_time)[0]
            prev_time = plt.date2num(
                plt.num2date(bry_time[0]) + relativedelta(days=-date_dt)
            )
            bry_time = np.append(prev_time, bry_time)
            del prev_time

        if nxt == 1 and len(bry_time) == 1:
            nxt_time = plt.date2num(plt.num2date(bry_time[-1]) + relativedelta(days=30))
            bry_time = np.append(bry_time, nxt_time)
            del nxt_time
        elif nxt == 1 and len(bry_time) > 1:
            date_dt = np.gradient(bry_time)[-1]
            nxt_time = plt.date2num(
                plt.num2date(bry_time[-1]) + relativedelta(days=date_dt)
            )
            bry_time = np.append(bry_time, nxt_time)
            del nxt_time

        nc = netcdf.Dataset(bdy_filename, "a")

        nc.Input_data_type = bry_def["inputdata"]
        nc.variables["bry_time"].cycle = bry_def["cycle_bry"]
        nc.variables["bry_time"][:] = bry_time
        if bry_def["cycle_bry"] == 0:
            nc.variables["bry_time"].units = (
                f"days since {origindate:%Y-%m-%d %H:%M:%S}"
            )
        # --- Loop on boundaries ------------------------------------------

        if len(bry_def["tracers"]) == 0:
            var_loop = ["ssh", "velocity"]
        else:
            var_loop = ["ssh", "tracers", "velocity"]

        for boundary, is_open in zip(
            bry_def["obc_dict"].keys(), bry_def["obc_dict"].values()
        ):
            if is_open:
                for lvars in var_loop:
                    print(
                        "\n     Processing *%s* for %sern boundary" % (lvars, boundary)
                    )
                    print("     ------------------------------------------")
                    if lvars == "ssh":
                        (zeta, NzGood) = interp_tools.interp_tracers(
                            inpdat,
                            lvars,
                            -1,
                            crocogrd,
                            dtmin,
                            dtmax,
                            prev,
                            nxt,
                            boundary[0].upper(),
                        )
                        z_rho = crocogrd.scoord2z_r(zeta=zeta, bdy="_" + boundary)
                        z_w = crocogrd.scoord2z_w(zeta=zeta, bdy="_" + boundary)

                    elif lvars == "tracers":
                        trac_dict = {}
                        for trc in bry_def["tracers"]:
                            print(f"\nIn tracers processing {trc}")
                            trac_dict[trc] = interp_tools.interp(
                                inpdat,
                                trc,
                                bry_def["Nzgoodmin"],
                                z_rho,
                                crocogrd,
                                dtmin,
                                dtmax,
                                prev,
                                nxt,
                                bdy=boundary[0].upper(),
                            )

                    elif lvars == "velocity":

                        cosa = np.cos(getattr(crocogrd, "angle_" + boundary))
                        sina = np.sin(getattr(crocogrd, "angle_" + boundary))

                        [
                            ubar_ogcm,
                            vbar_ogcm,
                            ubar,
                            vbar,
                        ] = interp_tools.compute_uvbar_ogcm(
                            inpdat,
                            cosa,
                            sina,
                            crocogrd,
                            dtmin,
                            dtmax,
                            prev,
                            nxt,
                            bdy=boundary[0].upper(),
                        )

                        [u, v] = interp_tools.interp_uv(
                            inpdat,
                            bry_def["Nzgoodmin"],
                            z_rho,
                            cosa,
                            sina,
                            crocogrd,
                            dtmin,
                            dtmax,
                            prev,
                            nxt,
                            bdy=boundary[0].upper(),
                        )

                        if bry_def["conserv"] == 1:
                            ubar_croco = sig_tools.vintegr4D(
                                u,
                                grd_tools.rho2u(z_w),
                                grd_tools.rho2u(z_rho),
                                np.nan,
                                np.nan,
                            )[0] / grd_tools.rho2u(getattr(crocogrd, "h_" + boundary))
                            vbar_croco = sig_tools.vintegr4D(
                                v,
                                grd_tools.rho2v(z_w),
                                grd_tools.rho2v(z_rho),
                                np.nan,
                                np.nan,
                            )[0] / grd_tools.rho2v(getattr(crocogrd, "h_" + boundary))

                            u = u - np.tile(
                                ubar_croco[:, np.newaxis, :, :],
                                (1, z_rho.shape[1], 1, 1),
                            )
                            u = u + np.tile(
                                ubar[:, np.newaxis, :, :], (1, z_rho.shape[1], 1, 1)
                            )

                            v = v - np.tile(
                                vbar_croco[:, np.newaxis, :, :],
                                (1, z_rho.shape[1], 1, 1),
                            )
                            v = v + np.tile(
                                vbar[:, np.newaxis, :, :], (1, z_rho.shape[1], 1, 1)
                            )

                # --- Saving in netcdf ------------------------------------------------

                print(f"\nSaving {boundary}ern boundary in Netcdf")
                print("----------------------------------")

                # handle indices (as 2 points where taken next to bdy)
                if str(boundary) == "west" and is_open:
                    indices3D = (
                        slice(None, None),
                        slice(None, None),
                        slice(None, None),
                        0,
                    )  # T,N,J,i=0
                    indices2D = (slice(None, None), slice(None, None), 0)  # T,J,i=0
                elif str(boundary) == "east" and is_open:
                    indices3D = (
                        slice(None, None),
                        slice(None, None),
                        slice(None, None),
                        -1,
                    )  # T,N,J,i=last
                    indices2D = (slice(None, None), slice(None, None), -1)  # T,J,i=last
                elif str(boundary) == "south" and is_open:
                    indices3D = (
                        slice(None, None),
                        slice(None, None),
                        0,
                        slice(None, None),
                    )  # T,N,j=0,I
                    indices2D = (slice(None, None), 0, slice(None, None))  # T,j=0,I
                elif str(boundary) == "north" and is_open:
                    indices3D = (
                        slice(None, None),
                        slice(None, None),
                        -1,
                        slice(None, None),
                    )  # T,N,j=last,I
                    indices2D = (slice(None, None), -1, slice(None, None))  # T,j=last,I

                mask_zet = np.tile(
                    getattr(crocogrd, "maskr_" + boundary), [zeta.shape[0], 1, 1]
                )
                if "velocity" in var_loop:
                    mask_u = np.tile(
                        getattr(crocogrd, "umask_" + boundary),
                        [u.shape[0], u.shape[1], 1, 1],
                    )
                    mask_v = np.tile(
                        getattr(crocogrd, "vmask_" + boundary),
                        [u.shape[0], u.shape[1], 1, 1],
                    )
                    mask_ubar = np.tile(
                        getattr(crocogrd, "umask_" + boundary), [u.shape[0], 1, 1]
                    )
                    mask_vbar = np.tile(
                        getattr(crocogrd, "vmask_" + boundary), [v.shape[0], 1, 1]
                    )

                nc.variables["zeta_" + str(boundary)][:] = (
                    zeta[indices2D] * mask_zet[indices2D]
                )
                nc.variables["u_" + str(boundary)][:] = u[indices3D] * mask_u[indices3D]
                nc.variables["v_" + str(boundary)][:] = v[indices3D] * mask_v[indices3D]
                nc.variables["ubar_" + str(boundary)][:] = (
                    ubar[indices2D] * mask_ubar[indices2D]
                )
                nc.variables["vbar_" + str(boundary)][:] = (
                    vbar[indices2D] * mask_vbar[indices2D]
                )

                if "tracers" in var_loop:
                    for varname, value in zip(trac_dict.keys(), trac_dict.values()):
                        mask_tra = np.tile(
                            getattr(crocogrd, "maskr_" + boundary),
                            [value.shape[0], value.shape[1], 1, 1],
                        )
                        nc.variables[f"{varname}_{boundary}"][:] = (
                            value[indices3D] * mask_tra[indices3D]
                        )

                # handle prev and nxt + save

        nc.close()
        # --- END writting netcdf ---------------------------------------------

        # --- Preparing time for next loop ------------------------------------
        # startloc = endloc + relativedelta(days=1, hours=-12)
        startloc = endloc
        if output_file_format == "DAILY":
            endloc = startloc + pandas.Timedelta(1, unit="day")
        elif output_file_format == "MONTHLY":
            endloc = startloc + relativedelta(months=1, days=-1, hour=12)
        elif output_file_format == "YEARLY":
            yearloc = plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S")
            if yearloc.year == int(bry_def["Yend"]):
                endloc = plt.num2date(dtend).replace(tzinfo=None)
            else:
                endloc = plt.datetime.datetime(int(yearloc.year), 12, 31, 12)
        elif bry_def["output_file_format"].upper() == "FULL":
            endloc = startloc

    return 0


if __name__ == "__main__":
    sys.exit(run_make_bry())
