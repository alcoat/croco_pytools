#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 15:30:02 2025

@author: annelou
"""

class aforc_class:
    def __init__(self):
        self.raw_name = {}
        self.units = {}
        self.conv_cff = {}
        self.iscumul =  {}
        self.filename = {}

    def add_var(self, key, value):
        self.raw_name[key] = value

    def set_conv(self, key, conv_cff):
        self.conv_cff[key] = conv_cff
        
    def set_cumul(self, key, iscumul):
        self.iscumul[key] = iscumul

    def file_name(self, key, value):
        self.filename[key] = value

    def get_var(self, key):
        return self.raw_name.get(key, None)

    def get_conv_cff(self, key):
        return self.conv_cff.get(key,None)
    
    def get_iscumul(self, key):
        return self.iscumul.get(key,None)

    def get_filename(self, key):
        return self.filename.get(key, None)

def create_class(var_info,multifiles):
    variables = aforc_class()
    for i in range(len(var_info)):
        variables.add_var(var_info[i][0],var_info[i][1])
        variables.set_conv(var_info[i][0],var_info[i][2])
        variables.set_cumul(var_info[i][0],var_info[i][3])
        if multifiles:
            variables.file_name(var_info[i][0],var_info[i][4])
    return variables
