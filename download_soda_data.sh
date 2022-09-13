#!/bin/bash
set -e
###
# The aim of this script is to download raw netcdf for SODA3.12.3 data.
# Here we just follow the instruction given here:
# https://www2.atmos.umd.edu/~ocean/index_files/soda3.12.2_mn_download.htm
# 
# If you want other SODA version you can have a look here:
# http://www.soda.umd.edu/
#
# Data are downloaded by years for montly data of every 5 days for 5-daily
# 1 year of monthly data is about 3.3-Go
# 5-days data are about 286-Mo
#
# Those data goes from 1980 to 2017
### 
YEAR_START=2017
MONTH_START=1

YEAR_END=2017
MONTH_END=1

#
### KIND OF DATA ###
kdata="MONTHLY" # 5DAILY or MONTHLY
####################
OUTDIR="./" # where to put data

###################################################################
########## END USER CHANGES #######################################
###################################################################

if [[ ${kdata} == "5DAILY" ]]; then
    prefix='soda3.12.2_5dy_ocean_reg_'
elif [[ ${kdata} == "MONTHLY" ]]; then
    prefix='soda3.12.2_mn_ocean_reg_'
else
    echo "Please specify what kind of data you want (5DAILY or MONTHLY), exit...."; exit 1
fi
###

### loop ###
for YEAR in `seq ${YEAR_START} ${YEAR_END}`; do
	if [[ ${kdata} == "MONTHLY" ]]; then
            wget -r -l1 --no-parent --progress=bar -nd -A.nc https://dsrs.atmos.umd.edu/DATA/soda3.12.2/REGRIDED/ocean/${prefix}${YEAR}.nc
			mv ${prefix}${YEAR}.nc ${OUTDIR}
	elif [[ ${kdata} == "5DAILY" ]]; then
        [[ ${YEAR} == ${YEAR_START} ]] && mstart=${MONTH_START} || mstart=1 
        [[ ${YEAR} == ${YEAR_END} ]] && mend=${MONTH_END} || mend=12
        for MONTH in `seq ${mstart} ${mend}`; do
	   	    wget -r -l1 --no-parent --progress=bar -nd -A.nc https://dsrs.atmos.umd.edu/DATA/soda3.12.2/REGRIDED/ocean/${prefix}${YEAR}_$(printf "%02d" ${MONTH})_*.nc
			 mv ${prefix}${YEAR}_$(printf "%02d" ${MONTH})_*.nc ${OUTDIR}
			if [[${MONTH} == ${MONTH_END} && ${YEAR} == ${YEAR_END}]]; then
			    echo "Download next month as it is not monthly separated"		
				wget -r -l1 --no-parent --progress=bar -nd -A.nc https://dsrs.atmos.umd.edu/DATA/soda3.12.2/REGRIDED/ocean/${prefix}${YEAR}_$(printf "%02d" $(( ${MONTH}+1 )) )_*.nc
				mv ${prefix}${YEAR}_$(printf "%02d" $(( ${MONTH}+1 )))_*.nc ${OUTDIR}
            fi
        done
    fi  
done 
