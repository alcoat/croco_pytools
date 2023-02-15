# -*- coding: utf-8 -*-
#!/usr/bin/python

import numpy as np
import datetime as dt
import pandas as pd
import os
import shutil
import fileinput
import subprocess
import wget


class oper_carpetas():
    def __init__ (self,ruta_raiz,fecha_ini):
        self.ruta_runt = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/'
        self.ruta_data = f'{ruta_raiz}data/{fecha_ini.strftime("%Y%m%d%H")}/'
        self.ruta_wrf = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/'
        self.ruta_cpl = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/cpl/'
    def crear_carpetas_l1(self):
        subprocess.call(['mkdir',self.ruta_runt])
        subprocess.call(['mkdir',self.ruta_data])
    def crear_carpetas_l2(self):
        subprocess.call(['mkdir',f'{self.ruta_runt}wrf/'])
        subprocess.call(['mkdir',f'{self.ruta_runt}cpl/'])
    def crear_carpetas_l3(self):
        subprocess.call(['mkdir',f'{self.ruta_wrf}gfs/'])
        subprocess.call(['mkdir',f'{self.ruta_wrf}wps/'])
        subprocess.call(['mkdir',f'{self.ruta_wrf}arw/'])
        subprocess.call(['mkdir',f'{self.ruta_data}wrf/'])
        subprocess.call(['mkdir',f'{self.ruta_data}croco/'])
        subprocess.call(['mkdir',f'{self.ruta_data}cpl/'])
        subprocess.call(['mkdir',f'{self.ruta_data}gfs/'])
        subprocess.call(['mkdir',f'{self.ruta_cpl}croco_files/'])

class oper_verifica():
    def __init__ (self,ruta_raiz,fecha_ini,nh):
        self.ruta_data = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/gfs/'
        self.nh = nh+3
        self.fini = fecha_ini
    def verifica_gfs(self):
        file_list = list(map(lambda x: 'gfs.t{0}z.pgrb2.0p25.f{1:03}'.format(self.fini.strftime('%H'),x),range(0,self.nh,3)))
        archivos = os.listdir(f'{self.ruta_data}')
        check_list = []
        for gfs_archivo in file_list:
            ruta_archivo = os.path.join(self.ruta_data, gfs_archivo)
            print (gfs_archivo,os.path.exists(ruta_archivo))
            if os.path.exists(ruta_archivo) == True:
                check_list.append(ruta_archivo)
        return check_list

class oper_copiar():
    def __init__ (self,ruta_raiz,fecha_ini):
        self.ruta_gfs = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/gfs/'
        self.ruta_arw = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/arw/'
        self.ruta_data_wrf = f'{ruta_raiz}data/{fecha_ini.strftime("%Y%m%d%H")}/wrf/'
        self.ruta_data_gfs = f'{ruta_raiz}data/{fecha_ini.strftime("%Y%m%d%H")}/gfs/'
    def copia_wrfout(self):
        archivos = os.listdir(f'{self.ruta_arw}')
        for archivo in archivos:
            ruta_archivo = os.path.join(self.ruta_arw, archivo)
            if os.path.isfile(ruta_archivo):
                if ruta_archivo.split('/')[-1][0:6] == 'wrfout':
                    shutil.copy(ruta_archivo, self.ruta_data_wrf)
    def copia_wrfhr(self):
        archivos = os.listdir(f'{self.ruta_arw}')
        for archivo in archivos:
            ruta_archivo = os.path.join(self.ruta_arw, archivo)
            if os.path.isfile(ruta_archivo):
                if ruta_archivo.split('/')[-1][0:7] == 'wrfxtrm':
                    shutil.copy(ruta_archivo, self.ruta_data_wrf)
    def copia_wrfrst(self):
        archivos = os.listdir(f'{self.ruta_arw}')
        for archivo in archivos:
            ruta_archivo = os.path.join(self.ruta_arw, archivo)
            if os.path.isfile(ruta_archivo):
                if ruta_archivo.split('/')[-1][0:6] == 'wrfrst':
                    shutil.copy(ruta_archivo, self.ruta_data_wrf)
    def copia_config(self):
        archivos = os.listdir(f'{self.ruta_arw}')
        for archivo in archivos:
            ruta_archivo = os.path.join(self.ruta_arw, archivo)
            if os.path.isfile(ruta_archivo):
                if(ruta_archivo.split('/')[-1] == 'namelist.input') | (ruta_archivo.split('/')[-1] == 'wrfbdy_d01') | \
                  (ruta_archivo.split('/')[-1] == 'wrfinput_d01') | (ruta_archivo.split('/')[-1] == 'wrflowinp_d01'):
                    shutil.copy(ruta_archivo, self.ruta_data_wrf)
    def copia_gfs(self):
        archivos = os.listdir(f'{self.ruta_gfs}')
        for archivo in archivos:
            ruta_archivo = os.path.join(self.ruta_gfs, archivo)
            if os.path.isfile(ruta_archivo):
                if ruta_archivo.split('/')[-1][0:3] == 'gfs':
                    shutil.copy(ruta_archivo, self.ruta_data_gfs)

class wrf_gfs():
    def __init__ (self,ruta_raiz,fecha_ini,nh):
        self.ruta_gfs = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/gfs/'
        self.noaa_url = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/'
        self.date_str = fecha_ini.strftime('%Y%m%d/%H/atmos/')
        self.noaa_fld = f'gfs.{self.date_str}'
        self.fini = fecha_ini
        self.nh = nh+3
    def descarga_gfs(self):
        file_list = list(map(lambda x: 'gfs.t{0}z.pgrb2.0p25.f{1:03}'.format(self.fini.strftime('%H'),x),range(0,self.nh,3)))
        print('hola')
        for file_name in file_list:
            print (f'{self.noaa_url}{self.noaa_fld}{file_name}')
            wget.download(f'{self.noaa_url}{self.noaa_fld}{file_name}', out=self.ruta_gfs)

class wrf_geo():
    def __init__ (self,ruta_raiz,ruta_geox,fecha_ini):
        self.ruta_info = f'{ruta_raiz}info/wrf/wps/'
        self.ruta_name = f'{ruta_raiz}info/wrf/name/'
        self.ruta_geog = f'{ruta_raiz}info/wrf/wps/WPS_GEOG/'
        self.ruta_wps = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/wps/'
        self.ruta_geox = ruta_geox
        self.fini = fecha_ini
    def copiar_archivos(self):
        shutil.copy(f'{self.ruta_info}GEOGRID.TBL', self.ruta_wps)
        shutil.copy(f'{self.ruta_name}namelist.wps_code', f'{self.ruta_wps}namelist.wps')
    def editar_namelist(self):
        with fileinput.FileInput(f'{self.ruta_wps}namelist.wps', inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('RUTA_GEOG', self.ruta_geog),end='')
    def ejecutar_geogrid(self):
        subprocess.call([f'{self.ruta_geox}geogrid.exe'],cwd=self.ruta_wps)

class wrf_lng():
    def __init__ (self,ruta_raiz,ruta_lngx,fecha_ini,suffix):
        self.ruta_gfs = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/gfs/'
        self.ruta_wps = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/wps/'
        self.ruta_lngx = ruta_lngx
        self.suffix = suffix
    def ejecutar_linkgrib(self):
        subprocess.run([f'{self.ruta_lngx}link_grib.csh', f'{self.ruta_gfs}{self.suffix}*'],cwd=self.ruta_wps)
        
class wrf_ung():
    def __init__ (self,ruta_raiz,ruta_ungx,fecha_ini,fecha_fin):
        self.ruta_info = f'{ruta_raiz}info/wrf/wps/'
        self.ruta_name = f'{ruta_raiz}info/wrf/name/'
        self.ruta_wps = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/wps/'
        self.ruta_ungx = ruta_ungx
        self.fini = fecha_ini
        self.ffin = fecha_fin
    def copiar_archivos(self):
        shutil.copy(f'{self.ruta_name}namelist.wps_code', f'{self.ruta_wps}namelist.wps')
        shutil.copy(f'{self.ruta_info}Vtable.GFS', f'{self.ruta_wps}Vtable')
    def editar_namelist(self):
        with fileinput.FileInput(f'{self.ruta_wps}namelist.wps', inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('FECHA_INI', self.fini.strftime("%Y-%m-%d_%H:%M:%S")),end='')
        with fileinput.FileInput(f'{self.ruta_wps}namelist.wps', inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('FECHA_FIN', self.ffin.strftime("%Y-%m-%d_%H:%M:%S")),end='')
        with fileinput.FileInput(f'{self.ruta_wps}namelist.wps', inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('TIPO-UNG', 'GFS'),end='')
    def ejecutar_ungrib(self):
        subprocess.call([f'{self.ruta_ungx}ungrib.exe'],cwd=self.ruta_wps)
    
class wrf_met():
    def __init__ (self,ruta_raiz,ruta_metx,fecha_ini,fecha_fin):
        self.ruta_info = f'{ruta_raiz}info/wrf/wps/'
        self.ruta_name = f'{ruta_raiz}info/wrf/name/'
        self.ruta_wps = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/wps/'
        self.ruta_metx = ruta_metx
        self.fini = fecha_ini
        self.ffin = fecha_fin
    def copiar_archivos(self):
        shutil.copy(f'{self.ruta_info}METGRID.TBL', f'{self.ruta_wps}')
        shutil.copy(f'{self.ruta_name}namelist.wps_code', f'{self.ruta_wps}namelist.wps')
    def editar_namelist(self):
        with fileinput.FileInput(f'{self.ruta_wps}namelist.wps', inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('FECHA_INI', self.fini.strftime("%Y-%m-%d_%H:%M:%S")),end='')
        with fileinput.FileInput(f'{self.ruta_wps}namelist.wps', inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('FECHA_FIN', self.ffin.strftime("%Y-%m-%d_%H:%M:%S")),end='')
        with fileinput.FileInput(f'{self.ruta_wps}namelist.wps', inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('TIPO-MET', 'GFS'),end='')
    def ejecutar_metgrid(self):
        subprocess.call([f'{self.ruta_metx}metgrid.exe'],cwd=self.ruta_wps)


class wrf_real():
    def __init__ (self,ruta_raiz,ruta_realx,fecha_ini,fecha_fin):
        self.ruta_info = f'{ruta_raiz}info/wrf/arw/'
        self.ruta_name = f'{ruta_raiz}info/wrf/name/'
        self.ruta_wps = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/wps/'
        self.ruta_arw = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/arw/'
        self.ruta_realx = ruta_realx
        self.fini = fecha_ini
        self.ffin = fecha_fin
    def copiar_archivos(self):
        shutil.copy(f'{self.ruta_name}namelist.input_code', f'{self.ruta_arw}namelist.input')
        archivos = os.listdir(f'{self.ruta_info}')
        for archivo in archivos:
            ruta_archivo = os.path.join(self.ruta_info, archivo)
            if os.path.isfile(ruta_archivo):
                shutil.copy(ruta_archivo, self.ruta_arw)

        archivos = os.listdir(f'{self.ruta_wps}')
        for archivo in archivos:
            ruta_archivo = os.path.join(self.ruta_wps, archivo)
            if os.path.isfile(ruta_archivo):
                if archivo.split('.')[0] == 'met_em':
                    shutil.copy(ruta_archivo, self.ruta_arw)
    def editar_namelist(self):
        with fileinput.FileInput(f'{self.ruta_arw}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('YI',self.fini.strftime("%Y")),end='')
        with fileinput.FileInput(f'{self.ruta_arw}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('MI',self.fini.strftime("%m")),end='')
        with fileinput.FileInput(f'{self.ruta_arw}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('DI',self.fini.strftime("%d")),end='')              
        with fileinput.FileInput(f'{self.ruta_arw}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('HI',self.fini.strftime("%H")),end='')
        with fileinput.FileInput(f'{self.ruta_arw}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('YF',self.ffin.strftime("%Y")),end='')
        with fileinput.FileInput(f'{self.ruta_arw}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('MF',self.ffin.strftime("%m")),end='')
        with fileinput.FileInput(f'{self.ruta_arw}namelist.input', inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('DF',self.ffin.strftime("%d")),end='')              
        with fileinput.FileInput(f'{self.ruta_arw}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('HF',self.ffin.strftime("%H")),end='')
        with fileinput.FileInput(f'{self.ruta_arw}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('HORA_PRONOSTICO,',str(int((self.ffin-self.fini).total_seconds()/3600.))),end='')
    def ejecutar_real(self):
        subprocess.call([f'{self.ruta_realx}real.exe'],cwd=self.ruta_arw)


class wrf_wrf():
    def __init__ (self,ruta_raiz,ruta_wrfx,ruta_mpix,nproc,name_nodo,fecha_ini):
        self.ruta_info = f'{ruta_raiz}info/wrf/arw/'
        self.ruta_arw = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/arw/'
        self.ruta_wrfx = ruta_wrfx
        self.ruta_mpix = ruta_mpix
        self.nproc = nproc
        self.name_nodo = name_nodo
    def ejecutar_wrf(self):
        subprocess.run([f'{self.ruta_mpix}','-np',f'{self.nproc}','-nameserver',
                        f'{self.name_nodo}',f'{self.ruta_wrfx}wrf.exe'],
                        cwd=self.ruta_arw)


