#!/usr/bin/env python
"""
===========================================================================
Further Information:
  http://www.croco-ocean.org

This file is part of CROCOTOOLS

Create a CROCO forcing file for tides
In the current state the script can handle:
    - TPXO (raw data)
    - FES2014

To add a new dataset you just have to go in Readers/tides_readers.py and
create a dico. Follow information given inside the file to fill it properly

The periodicity of a dataset (if it is -180/180 or 0/360) is handle by the
script. Just download the dataset and it should be fine

The script works as follow:
    - reads croco_grd
    - reads input data and restricts the area of spatial coverage
    - creates croco_frc.nc
    - Loop on choosen waves
    - Writes data in netcdf

===========================================================================
"""

__author__ = "Mathieu Le Corre"
__email__ = "mathieu.le.corre@shom.fr"
__date__ = "2022-11"
__license__ = "GPL3"


import argparse
import pathlib
import sys
from os import path
import yaml

import cftime
import netCDF4 as netcdf
import numpy as np

sys.path.append(pathlib.Path(__file__).resolve().parent / "Modules")
sys.path.append(pathlib.Path(__file__).resolve().parent / "Readers" )
import croco_class as Croco
import interp_tools
import tides_class as Inp


def get_args():
    """get arguments"""
    parser = argparse.ArgumentParser(
        description="make tides",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--config_file",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parent / "make_tides_def.yml",
        help="Config file",
    )
    args = parser.parse_args()

    return args


def run_make_tides():
    """make tides program"""
    args = get_args()

    with open(args.config_file, "r", encoding="utf8") as infile:
        tide_def = yaml.safe_load(infile)

    print(tide_def)

    # --- START MAIN SCRIPT ----------------------------------------------------

    # read tides and periods
    tides_param = np.loadtxt(
        pathlib.Path(__file__).resolve().parent / "Modules" / "tides.txt", skiprows=3, comments="^", usecols=(0, 4), dtype=str
    )
    tides_names = tides_param[:, 0].tolist()
    tides_periods = tides_param[:, 1].tolist()
    tides_names = [x.lower() for x in tides_names]
    tides_names = np.array(tides_names)

    # Load croco_grd
    sigma_params = {"theta_s": 0, "theta_b": 0, "N": 1, "hc": 1}  # no needs as 2d
    croco_grd_path = pathlib.Path(tide_def["croco_dir"]) / tide_def["croco_grd"]
    crocogrd = Croco.CROCO_grd(croco_grd_path, sigma_params)

    # --- Load input (restricted to croco_grd) ----------------------------

    if tide_def["multi_files"]:
        if tide_def["waves_separated"]:
            input_file_ssh = []
            input_file_u = []
            input_file_v = []
            for inp in tide_def["tides"]:
                if path.isfile(
                    tide_def["input_dir"]
                    + tide_def["elev_file"].replace("<tides>", inp)
                ):
                    input_file_ssh += [
                        tide_def["input_dir"]
                        + tide_def["elev_file"].replace("<tides>", inp)
                    ]
                elif path.isfile(
                    tide_def["input_dir"]
                    + tide_def["elev_file"].replace("<tides>", inp.lower())
                ):
                    input_file_ssh += [
                        tide_def["input_dir"]
                        + tide_def["elev_file"].replace("<tides>", inp.lower())
                    ]
                else:
                    el_fname = tide_def["input_dir"] + tide_def["elev_file"].replace(
                        "<tides>", inp
                    )
                    print(f"Elevation file {el_fname} for wave {inp} is missing")
                    sys.exit(1)

                if tide_def["cur"]:
                    if path.isfile(
                        tide_def["input_dir"]
                        + tide_def["u_file"].replace("<tides>", inp)
                    ):
                        input_file_u += [
                            tide_def["input_dir"]
                            + tide_def["u_file"].replace("<tides>", inp)
                        ]
                    elif path.isfile(
                        tide_def["input_dir"]
                        + tide_def["u_file"].replace("<tides>", inp.lower())
                    ):
                        input_file_u += [
                            tide_def["input_dir"]
                            + tide_def["u_file"].replace("<tides>", inp.lower())
                        ]
                    else:
                        print(f"Eastward current file for wave {inp} is missing")
                        sys.exit(1)

                    if path.isfile(
                        tide_def["input_dir"]
                        + tide_def["v_file"].replace("<tides>", inp)
                    ):
                        input_file_v += [
                            tide_def["input_dir"]
                            + tide_def["v_file"].replace("<tides>", inp)
                        ]
                    elif path.isfile(
                        tide_def["input_dir"]
                        + tide_def["v_file"].replace("<tides>", inp.lower())
                    ):
                        input_file_v += [
                            tide_def["input_dir"]
                            + tide_def["v_file"].replace("<tides>", inp.lower())
                        ]
                    else:
                        print(f"Northward current file for wave {inp} is missing")
                        sys.exit(1)

            input_file_ssh = list(input_file_ssh)
            if tide_def["cur"]:
                input_file_u = list(input_file_u)
                input_file_v = list(input_file_v)
            else:
                input_file_u = None
                input_file_v = None
        else:
            input_file_ssh = list(tide_def["input_dir"] + tide_def["elev_file"])
            if tide_def["cur"]:
                input_file_u = list(tide_def["input_dir"] + tide_def["u_file"])
                input_file_v = list(tide_def["input_dir"] + tide_def["v_file"])
            else:
                input_file_u = None
                input_file_v = None
    else:
        input_file_ssh = list([tide_def["input_dir"] + tide_def["input_file"]])
        if tide_def["cur"]:
            input_file_u = list([tide_def["input_dir"] + tide_def["input_file"]])
            input_file_v = list([tide_def["input_dir"] + tide_def["input_file"]])
        else:
            input_file_u = None
            input_file_v = None

    inpdat = Inp.getdata(
        tide_def["inputdata"],
        input_file_ssh,
        crocogrd,
        tide_def["input_type"],
        tide_def["tides"],
        input_file_u,
        input_file_v,
    )

    # --- Create the initial file -----------------------------------------
    croco_tide_pname = pathlib.Path(tide_def["croco_dir"]) / tide_def["croco_filename"]
    print(f"Create file {croco_tide_pname}")
    Croco.CROCO.create_tide_nc(
        None,
        croco_tide_pname,
        crocogrd,
        cur=tide_def["cur"],
        pot=tide_def["pot"],
    )

    date = None
    if tide_def["Correction_ssh"] or tide_def["Correction_uv"]:
        date = cftime.datetime(tide_def["Yini"], tide_def["Mini"], tide_def["Dini"])
        date_orig = cftime.datetime(
            tide_def["Yorig"], tide_def["Morig"], tide_def["Dorig"]
        )

    todo = ["H"]

    if tide_def["cur"]:
        todo += ["cur"]
    coslat2 = None
    sin2lat = None
    if tide_def["pot"]:
        todo += ["pot"]
        coslat2 = np.cos(np.deg2rad(crocogrd.lat)) ** 2
        sin2lat = np.sin(2.0 * np.deg2rad(crocogrd.lat))

    # --- Start loop on waves --------------------------------------------
    print(f"Fill file %s", croco_tide_pname)
    nc = netcdf.Dataset(croco_tide_pname, "a")

    if tide_def["Correction_ssh"] or tide_def["Correction_uv"]:
        nc.Nodal_Correction = "".join(("Origin time is ", str(date_orig)))
    else:
        nc.Nodal_Correction = "No nodal correction"

    for i, tide in enumerate(tide_def["tides"]):
        print(f"\nProcessing *{tide}* wave")
        print("-----------------------")
        index = np.argwhere(tides_names == tide.lower())
        if len(index) > 0:
            print(f"  tides {tide} is in the list")
            # get period
            index = index[0][0]
            period = float(tides_periods[index])
            print(f"  Period of the wave {tide} is {period}")

            nc.variables["tide_period"][i] = period
        if tide_def["multi_files"]:
            # In this case waves had been concatenated in the order they appear in the list
            tndx = i
        else:
            # read ntime/periods dimension and find the closest wave
            tndx = np.argwhere(abs(inpdat.ntides - period) < 1e-3)
            if len(tndx) == 0:
                print(f"  Did not find wave {tide} in input file")
                sys.exit(1)
            elif len(tndx) > 1:  # increase the precision
                tndx = np.argwhere(abs(inpdat.ntides - period) < 1e-4)
                if len(tndx) == 0:
                    print("  Did not find wave {tide} in input file")
                    sys.exit(1)
                else:
                    tndx = tndx[0]
            else:
                tndx = tndx[0]

        [pf, pu, mk_b] = inpdat.egbert_correction(tide, date)
        # For phase shift time should be in seconds relatively Jan 1 1992
        # As mk_b is the wave phase at this date
        t0 = cftime.date2num(date_orig, "seconds since 1992-01-01:00:00:00")
        if tide_def["Correction_ssh"] or tide_def["Correction_uv"]:
            correc_amp = pf
            correc_phase = mk_b + np.deg2rad(t0 / (period * 10)) + pu
            #              |--- phase at origin time ---|  |nodal cor at ini time|
        else:
            correc_amp = 1
            correc_phase = mk_b + np.deg2rad(t0 / (period * 10))

        # --- Start loop on var ------------------------------------------
        for lvars in todo:
            # get data
            if lvars == "H":
                print("\n  Processing tidal elevation")
                print("  -------------------------")
                (tide_complex, nz_good) = interp_tools.interp_tides(
                    inpdat, lvars, -1, crocogrd, tndx, tndx, tide_def["input_type"]
                )
                if tide_def["Correction_ssh"]:
                    tide_amp = np.ma.abs(tide_complex) * correc_amp
                    if "tpxo" in tide_def["inputdata"]:
                        tide_phase = np.mod(
                            np.ma.angle(tide_complex) * -180 / np.pi
                            - correc_phase * 180 / np.pi,
                            360,
                        )
                    else:
                        tide_phase = np.mod(
                            np.ma.angle(tide_complex) * 180.0 / np.pi
                            - correc_phase * 180 / np.pi,
                            360,
                        )
                else:
                    tide_amp = np.ma.abs(tide_complex)
                    if "tpxo" in tide_def["inputdata"]:
                        tide_phase = np.mod(
                            np.ma.angle(tide_complex) * -180 / np.pi, 360
                        )
                    else:
                        tide_phase = np.mod(
                            np.ma.angle(tide_complex) * 180.0 / np.pi, 360
                        )

                nc.variables["tide_Ephase"][i, :] = tide_phase * crocogrd.maskr
                nc.variables["tide_Eamp"][i, :] = tide_amp * crocogrd.maskr

            #########################
            elif lvars == "cur":
                print("\n  Processing tidal currents")
                print("  -------------------------")
                (u_tide_complex, v_tide_complex, nz_good) = interp_tools.interp_tides(
                    inpdat, lvars, -1, crocogrd, tndx, tndx, tide_def["input_type"]
                )

                if tide_def["Correction_uv"]:
                    u_tide_amp = np.ma.abs(u_tide_complex) * correc_amp
                    v_tide_amp = np.ma.abs(v_tide_complex) * correc_amp

                    if "tpxo" in tide_def["inputdata"]:
                        u_tide_phase = np.mod(
                            np.ma.angle(u_tide_complex) * -180 / np.pi
                            - correc_phase * 180 / np.pi,
                            360,
                        )
                        v_tide_phase = np.mod(
                            np.ma.angle(v_tide_complex) * -180 / np.pi
                            - correc_phase * 180 / np.pi,
                            360,
                        )
                    else:
                        u_tide_phase = np.mod(
                            np.ma.angle(u_tide_complex) * 180.0 / np.pi
                            - correc_phase * 180 / np.pi,
                            360,
                        )
                        v_tide_phase = np.mod(
                            np.ma.angle(v_tide_complex) * 180.0 / np.pi
                            - correc_phase * 180 / np.pi,
                            360,
                        )
                else:
                    u_tide_amp = np.ma.abs(u_tide_complex)
                    v_tide_amp = np.ma.abs(v_tide_complex)

                    if "tpxo" in tide_def["inputdata"]:
                        u_tide_phase = np.mod(
                            np.ma.angle(u_tide_complex) * -180 / np.pi, 360
                        )
                        v_tide_phase = np.mod(
                            np.ma.angle(v_tide_complex) * -180 / np.pi, 360
                        )
                    else:
                        u_tide_phase = np.mod(
                            np.ma.angle(u_tide_complex) * 180.0 / np.pi, 360
                        )
                        v_tide_phase = np.mod(
                            np.ma.angle(v_tide_complex) * 180.0 / np.pi, 360
                        )

                major, eccentricity, inclination, phase = inpdat.ap2ep(
                    u_tide_amp, u_tide_phase, v_tide_amp, v_tide_phase
                )

                nc.variables["tide_Cmin"][i, :, :] = (
                    major[:, :] * eccentricity[:, :] * crocogrd.maskr
                )
                nc.variables["tide_Cmax"][i, :, :] = major[:, :] * crocogrd.maskr
                nc.variables["tide_Cangle"][i, :, :] = (
                    inclination[:, :] * crocogrd.maskr
                )
                nc.variables["tide_Cphase"][i, :, :] = phase[:, :] * crocogrd.maskr
            #########################
            elif lvars == "pot":
                print("\n  Processing equilibrium tidal potential")
                print("  --------------------------------------")
                try:
                    coef = getattr(inpdat.pot_tide, tide.lower())
                except:
                    try:
                        # some waves start with a number (ex: 2N2) and python do not like it
                        coef = getattr(inpdat.pot_tide, tide.lower())
                    except:
                        print("No potential prameter defined for wave {input_wav}")
                        coef = [1, 0]

                if period < 13:  # semidiurnal
                    p_amp = correc_amp * coef[0] * coef[1] * coslat2
                    p_pha = np.mod(-2 * crocogrd.lon - correc_phase * 180 / np.pi, 360)
                elif period < 26:  # diurnal
                    p_amp = correc_amp * coef[0] * coef[1] * sin2lat
                    p_pha = np.mod(-crocogrd.lon - correc_phase * 180 / np.pi, 360)
                else:  # long-term
                    p_amp = correc_amp * coef[0] * coef[1] * (1 - 1.5 * coslat2)
                    p_pha = np.mod(-correc_phase * 180 / np.pi, 360.0)

                nc.variables["tide_Pamp"][i, :, :] = p_amp * crocogrd.maskr
                nc.variables["tide_Pphase"][i, :, :] = p_pha * crocogrd.maskr

    nc.close()
    print("End processing")
    return 0


if __name__ == "__main__":
    sys.exit(run_make_tides())
