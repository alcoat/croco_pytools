#!/bin/bash
#
# Andres Sepulveda (DGEO - UdeC) 16/06/2024
# hectorsepulveda@udec.cl
#

#
# Find a better way to handle time
# 10901 = 28/08/2022
#        Timesteps are every 3 hours
#
# Add +/- 2 degrees to GRD area coverage
#
for i in {10901..11651}
do 
   ncks -v water_temp,salinity,surf_el,water_v,water_u -d time,$i,$i -d lat,-40.1,-23.6 -d lon,4.0,26.0 http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0 hycom_$i.nc
   echo $i
   sleep 0.5
done

for i in {10901..11651}
do 
   ncks --mk_rec_dmn time hycom_TMC_$i.nc new_hycom_$i.nc
   echo $i
done

ncrcat new_hycom_?????.nc hycom_all.nc

#
# Conversion from 181-360 should be handled by ibc_class.py but it doesn't. Fix
#
ncap2 -O -s 'where(lon > 180) lon=lon-360'   hycom_all.nc  hycom_all.nc

#
# Add if all ok
#
#rm hycom_?????.nc
#rm new_hycom_?????.nc

