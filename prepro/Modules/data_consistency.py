#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 20:02:29 2024

@author: menkes
"""
import sys
import numpy as np
import pylab as plt
# -----------------------------------------------------------
# FUNCTIONS USED BY make_aforc.py TO CHECK DATA CONSISTENCY
# -----------------------------------------------------------
# data : xarray.DataArray of one variable
def consistency(data):
    if data.dims['time']==0 :
        print('\nData is missing - No Data for the subset')
        sys.exit()
    if np.nanvar(np.gradient(plt.date2num(data['time'].values))) >=5: # Abnormal distribution of days
        Question = input( "Abnormal distribution of days (variance to high) \
          \nThis may be due to the use of different temproral resolution dataset.\
               \n Do you want to proceed?: y,[n] ") or 'no'
        if Question.lower() == ("n") or Question.lower() == ("no"):
            print('Aborting')
            sys.exit()
                    
