__author__ = "Mathieu Le Corre"
__email__ = "mathieu.le.corre@shom.fr"
__date__ = "2022-11"
__license__ = "GPL3"

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

# --- Dependencies ---------------------------------------------------------

import glob
import sys
from os import path

import cftime
import netCDF4 as netcdf
import numpy as np

sys.path.append("./Modules/")
sys.path.append("./Readers/")
import croco_class as Croco
import interp_tools
import tides_class as Inp

# --- USER CHANGES ---------------------------------------------------------
tide_def = {
    # Dates
    # Initial date
    "Yini": 2013,
    "Mini": 1,
    "Dini": 1,
    # Origin year
    "Yorig": 2000,  # 1900,1,1 if TIDES_MAS OR ANA_INITIAL+USE_CALENDAR defined in cppdef.h
    "Morig": 1,  # 1900,1,1 if TIDES_MAS OR ANA_INITIAL+USE_CALENDAR defined in cppdef.h
    "Dorig": 1,  # 1900,1,1 if TIDES_MAS OR ANA_INITIAL+USE_CALENDAR defined in cppdef.h
    # Input data information and formating
    # Note: if you are using a tpxo dataset please be sure to have somewhere in
    #       inputdata 'tpxo'. This will allow the code to use the OTIS (TPXO is obtained with it)
    #       convention a-b*i for complex.
    #       Also, if you have already preprocess TPXO, make sure that you have correct units and
    #       u,v are in m/s and not m**2/s
    "inputdata": "tpxo7_croco",  # Input data dictionnary as defined in the Readers/tides_reader.py
    "input_dir": "/home/opsys/DATA/TIDE/OSU/",
    "input_file": "TPXO7.nc",  # Leave empty if you have multiple files
    "input_type": "Re_Im",  # Format of the input data 'Amp_Phase'or 'Re_Im'
    "multi_files": False,  # Set to True if several input files
    "waves_separated": False,  # Set to True if input files waves are separated
    "elev_file": "h_<tides>_tpxo9_atlas_30_v5.nc",  # elevation file names. if wave_separated put <tides> where wave name is found
    "u_file": "u_<tides>_tpxo9_atlas_30_v5.nc",  # eastward currents file names. if wave_separated put <tides> where wave name is found
    "v_file": "u_<tides>_tpxo9_atlas_30_v5.nc",  # northward currents file names. if wave_separated put <tides> where wave name is found
    # CROCO grid informations
    "croco_dir": "../",
    "croco_grd": "croco_grd.nc",
    # Tide file informations
    "croco_filename": "croco_frc.nc",
    "tides": ["M2", "S2", "N2", "K2", "K1", "O1", "P1", "Q1", "Mf", "Mm"],
    "cur": True,  # Set to True if you to compute currents
    "pot": True,  # Set to True if you to compute potiential tides
    # Nodal correction
    "Correction_ssh": True,
    "Correction_uv": True,
}

print(tide_def)
# --- END USER CHANGES -----------------------------------------------------

# --- START MAIN SCRIPT ----------------------------------------------------

# read tides and periods
tides_param = np.loadtxt(
    "./Modules/tides.txt", skiprows=3, comments="^", usecols=(0, 4), dtype=str
)
tides_names = tides_param[:, 0].tolist()
tides_periods = tides_param[:, 1].tolist()
tides_names = [x.lower() for x in tides_names]
tides_names = np.array(tides_names)

# Load croco_grd
sigma_params = dict(theta_s=0, theta_b=0, N=1, hc=1)  # no needs as 2d
crocogrd = Croco.CROCO_grd(
    "".join((tide_def["croco_dir"], tide_def["croco_grd"])), sigma_params
)

# --- Load input (restricted to croco_grd) ----------------------------

if tide_def["multi_files"]:
    if tide_def["waves_separated"]:
        input_file_ssh = []
        input_file_u = []
        input_file_v = []
        for inp in tide_def['tides']:
            if path.isfile(
                tide_def["input_dir"] + tide_def["elev_file"].replace("<tides>", inp)
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
                sys.exit(
                    "Elevation file %s for wave %s is missing"
                    % (
                        tide_def["input_dir"]
                        + tide_def["elev_file"].replace("<tides>", inp),
                        inp,
                    )
                )

            if tide_def['cur']:
                if path.isfile(
                    tide_def["input_dir"] + tide_def["u_file"].replace("<tides>", inp)
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
                    sys.exit("Eastward current file for wave %s is missing" % inp)

                if path.isfile(
                    tide_def["input_dir"] + tide_def["v_file"].replace("<tides>", inp)
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
                    sys.exit("Northward current file for wave %s is missing" % inp)

        input_file_ssh = list(input_file_ssh)
        if tide_def['cur']:
            input_file_u = list(input_file_u)
            input_file_v = list(input_file_v)
        else:
            input_file_u = None
            input_file_v = None
    else:
        input_file_ssh = list(tide_def["input_dir"] + tide_def["elev_file"])
        if tide_def['cur']:
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

Croco.CROCO.create_tide_nc(
    None,
    "".join((tide_def["croco_dir"] + tide_def["croco_filename"])),
    crocogrd,
    cur=tide_def["cur"],
    pot=tide_def["pot"],
)

if tide_def["Correction_ssh"] or tide_def["Correction_uv"]:
    date = cftime.datetime(tide_def["Yini"], tide_def["Mini"], tide_def["Dini"])
    date_orig = cftime.datetime(tide_def["Yorig"], tide_def["Morig"], tide_def["Dorig"])

todo = ["H"]

if tide_def["cur"]:
    todo += ["cur"]
if tide_def["pot"]:
    todo += ["pot"]
    coslat2 = np.cos(np.deg2rad(crocogrd.lat)) ** 2
    sin2lat = np.sin(2.0 * np.deg2rad(crocogrd.lat))

# --- Start loop on waves --------------------------------------------
nc = netcdf.Dataset(tide_def["croco_dir"] + tide_def["croco_filename"], "a")

if tide_def["Correction_ssh"] or tide_def["Correction_uv"]:
    nc.Nodal_Correction = "".join(("Origin time is ", str(date_orig)))
else:
    nc.Nodal_Correction = "No nodal correction"

for i, tide in enumerate(tide_def["tides"]):
    print("\nProcessing *%s* wave" % (tide))
    print("-----------------------")
    index = np.argwhere(tides_names == tide.lower())
    if len(index) > 0:
        print("  tides %s is in the list" % (tide))
        # get period
        index = index[0][0]
        period = float(tides_periods[index])
        print("  Period of the wave %s is %f" % (tide, period))

        nc.variables["tide_period"][i] = period
    if tide_def["multi_files"]:
        # In this case waves had been concatenated in the order they appear in the list
        tndx = i
    else:
        # read ntime/periods dimension and find the closest wave
        tndx = np.argwhere(abs(inpdat.ntides - period) < 1e-3)
        if len(tndx) == 0:
            sys.exit("  Did not find wave %s in input file" % tide)
        elif len(tndx) > 1:  # increase the precision
            tndx = np.argwhere(abs(inpdat.ntides - period) < 1e-4)
            if len(tndx) == 0:
                sys.exit("  Did not find wave %s in input file" % tide)
            else:
                tndx = tndx[0]
        else:
            tndx = tndx[0]

    [pf, pu, mkB] = inpdat.egbert_correction(tide, date)
    # For phase shift time should be in seconds relatively Jan 1 1992
    # As mkB is the wave phase at this date
    t0 = cftime.date2num(date_orig, "seconds since 1992-01-01:00:00:00")
    if tide_def["Correction_ssh"] or tide_def["Correction_uv"]:
        correc_amp = pf
        correc_phase = mkB + np.deg2rad(t0 / (period * 10)) + pu
        #              |--- phase at origin time ---|  |nodal cor at ini time|
    else:
        correc_amp = 1
        correc_phase = mkB + np.deg2rad(t0 / (period * 10))

    # --- Start loop on var ------------------------------------------
    for lvars in todo:
        # get data
        if lvars == "H":
            print("\n  Processing tidal elevation")
            print("  -------------------------")
            (tide_complex, NzGood) = interp_tools.interp_tides(
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
                    tide_phase = np.mod(np.ma.angle(tide_complex) * -180 / np.pi, 360)
                else:
                    tide_phase = np.mod(np.ma.angle(tide_complex) * 180.0 / np.pi, 360)

            nc.variables["tide_Ephase"][i, :] = tide_phase * crocogrd.maskr
            nc.variables["tide_Eamp"][i, :] = tide_amp * crocogrd.maskr

        #########################
        elif lvars == "cur":
            print("\n  Processing tidal currents")
            print("  -------------------------")
            (u_tide_complex, v_tide_complex, NzGood) = interp_tools.interp_tides(
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
            nc.variables["tide_Cangle"][i, :, :] = inclination[:, :] * crocogrd.maskr
            nc.variables["tide_Cphase"][i, :, :] = phase[:, :] * crocogrd.maskr
        #########################
        elif lvars == "pot":
            print("\n  Processing equilibrium tidal potential")
            print("  --------------------------------------")
            try:
                coef = eval("".join(("inpdat.pot_tide.", tide.lower())))
            except:
                try:
                    # some waves start with a number (ex: 2N2) and python do not like it
                    coef = eval("".join(("inpdat.pot_tide._", tide.lower())))
                except:
                    print("No potential prameter defined for wave %s" % input_wav)
                    coef = [1, 0]

            if period < 13:  # semidiurnal
                Pamp = correc_amp * coef[0] * coef[1] * coslat2
                Ppha = np.mod(-2 * crocogrd.lon - correc_phase * 180 / np.pi, 360)
            elif period < 26:  # diurnal
                Pamp = correc_amp * coef[0] * coef[1] * sin2lat
                Ppha = np.mod(-crocogrd.lon - correc_phase * 180 / np.pi, 360)
            else:  # long-term
                Pamp = correc_amp * coef[0] * coef[1] * (1 - 1.5 * coslat2)
                Ppha = np.mod(-correc_phase * 180 / np.pi, 360.0)

            nc.variables["tide_Pamp"][i, :, :] = Pamp * crocogrd.maskr
            nc.variables["tide_Pphase"][i, :, :] = Ppha * crocogrd.maskr

nc.close()
