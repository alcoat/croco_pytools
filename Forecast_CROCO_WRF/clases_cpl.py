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
import netCDF4
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.interpolate import griddata
import pygrib


class oper_copiar():
    def __init__ (self,ruta_raiz,fecha_ini):
        self.ruta_info = f'{ruta_raiz}info/cpl/'
        self.ruta_runt = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/'
        self.ruta_data = f'{ruta_raiz}data/{fecha_ini.strftime("%Y%m%d%H")}/'
        self.ruta_runt_cpl = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/cpl/'
    def copia_wrf(self):
        archivos = os.listdir(f'{self.ruta_info}wrf/')
        for archivo in archivos:
            ruta_archivo = os.path.join(f'{self.ruta_info}wrf/', archivo)
            if os.path.isfile(ruta_archivo):
                    shutil.copy(ruta_archivo, self.ruta_runt_cpl)
    def copia_croco(self):
        subprocess.call(['mkdir',f'{self.ruta_runt_cpl}croco_files/'])
        subprocess.call(['mkdir',f'{self.ruta_runt_cpl}croco_outfl/'])
        archivos = os.listdir(f'{self.ruta_info}croco/')
        for archivo in archivos:
            ruta_archivo = os.path.join(f'{self.ruta_info}croco/', archivo)
            if os.path.isfile(ruta_archivo):
                if ruta_archivo.split('/')[-1][6:9] == 'grd':
                    shutil.copy(ruta_archivo, f'{self.ruta_runt_cpl}croco_files/')
                else:
                    shutil.copy(ruta_archivo, self.ruta_runt_cpl)
    def copia_oasis(self):
        archivos = os.listdir(f'{self.ruta_info}oasis/')
        for archivo in archivos:
            ruta_archivo = os.path.join(f'{self.ruta_info}oasis/', archivo)
            if os.path.isfile(ruta_archivo):
                    shutil.copy(ruta_archivo, self.ruta_runt_cpl)

class pre_wrf():
    def __init__ (self,ruta_raiz,fecha_ini,fecha_fin):
        self.ruta_runt_wrf = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/wrf/arw/'
        self.ruta_runt_cpl = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/cpl/'
        self.fini = fecha_ini
        self.ffin = fecha_fin
    def copia_wrf_archivos(self):
        archivos = os.listdir(f'{self.ruta_runt_wrf}')
        for archivo in archivos:
            ruta_archivo = os.path.join(self.ruta_runt_wrf, archivo)
            if os.path.isfile(ruta_archivo):
                if(ruta_archivo.split('/')[-1] == 'wrfbdy_d01') | \
                  (ruta_archivo.split('/')[-1] == 'wrfinput_d01') | \
                  (ruta_archivo.split('/')[-1] == 'wrflowinp_d01'):
                    shutil.copy(ruta_archivo, self.ruta_runt_cpl)
    def editar_namelist(self):
        with fileinput.FileInput(f'{self.ruta_runt_cpl}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('YI',self.fini.strftime("%Y")),end='')
        with fileinput.FileInput(f'{self.ruta_runt_cpl}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('MI',self.fini.strftime("%m")),end='')
        with fileinput.FileInput(f'{self.ruta_runt_cpl}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('DI',self.fini.strftime("%d")),end='')              
        with fileinput.FileInput(f'{self.ruta_runt_cpl}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('HI',self.fini.strftime("%H")),end='')
        with fileinput.FileInput(f'{self.ruta_runt_cpl}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('YF',self.ffin.strftime("%Y")),end='')
        with fileinput.FileInput(f'{self.ruta_runt_cpl}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('MF',self.ffin.strftime("%m")),end='')
        with fileinput.FileInput(f'{self.ruta_runt_cpl}namelist.input', inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('DF',self.ffin.strftime("%d")),end='')              
        with fileinput.FileInput(f'{self.ruta_runt_cpl}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('HF',self.ffin.strftime("%H")),end='')
        with fileinput.FileInput(f'{self.ruta_runt_cpl}namelist.input',inplace=True, backup='.bak') as file:
            for line in file:
                print (line.replace('HORA_PRONOSTICO,',str(int((self.ffin-self.fini).total_seconds()/3600.))),end='')

    def modifica_wrfinput(self):
        name_grd = f'{self.ruta_runt_cpl}/croco_files/croco_grd.nc'
        name_win = f'{self.ruta_runt_cpl}/wrfinput_d01'

        nc_croco = netCDF4.Dataset(name_grd,'r')
        mask_croco = nc_croco.variables['mask_rho'][:,:]
        lat_croco = nc_croco.variables['lat_rho'][:,:]
        lon_croco = nc_croco.variables['lon_rho'][:,:]
        nc_croco.close()

        nc_wrfi = netCDF4.Dataset(name_win,'r')
        lat_wrf = nc_wrfi.variables['XLAT'][0,:,:]
        lon_wrf = nc_wrfi.variables['XLONG'][0,:,:]
        nc_wrfi.close()

        valuesm = np.reshape(mask_croco,-1)
        pointsm = np.zeros((mask_croco.shape[0]*mask_croco.shape[1],2))
        xx = np.reshape(lon_croco,-1); yy = np.reshape(lat_croco,-1)
        pointsm[:,0] = xx; pointsm[:,1] = yy
        mask_croco_int = griddata(pointsm, valuesm, (lon_wrf, lat_wrf), method='nearest')

        nc_wrfi = netCDF4.Dataset(name_win,'a')
        nc_wrfi.variables['CPLMASK'][0,0,:,:] = mask_croco_int
        nc_wrfi.close()

    def grafica_mascara(self):
        name_grd = f'{self.ruta_runt_cpl}/croco_files/croco_grd.nc'
        name_win = f'{self.ruta_runt_cpl}/wrfinput_d01'
        nc_wrfi = netCDF4.Dataset(name_win,'r')
        cpl = nc_wrfi.variables['CPLMASK'][:,:,:,:]
        lat_wrf = nc_wrfi.variables['XLAT'][0,:,:]
        lon_wrf = nc_wrfi.variables['XLONG'][0,:,:]
        nc_wrfi.close()
        cmap = colors.ListedColormap(['white', 'darkblue'])
        bounds=[0,.5,1]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        fig = plt.figure(figsize=(10,6))
        ax = fig.add_subplot(111)
        mapa = ax.contourf(lon_wrf,lat_wrf,cpl[0,0,:,:],cmap=cmap, norm=norm, boundaries=bounds)
        plt.text(261.3-360,33.2,f'WRFINPUT-Mask Modified {self.fini.strftime("%Y %b %d %H:00")}',fontsize=20,fontweight='bold')
        cbar_ax = fig.add_axes([0.92, 0.14, 0.03, 0.68])
        cbar1 = plt.colorbar(mapa,cax=cbar_ax,orientation='vertical',ticks=[0,1])
        cbar1.set_label(r'Mask $[Land = 0, Ocean = 1]$',fontsize=18)
        cbar1.ax.tick_params(labelsize=14)
        plt.savefig(f'{self.ruta_runt_cpl}mask_mod.png',dpi=200,bbox_inches="tight")
        
class pre_oasis():
    def __init__ (self,ruta_raiz,fecha_ini,ruta_croco_his,base_croco,ruta_wrst,ruta_gfs,gfs_file):
        self.ruta_pros = f'{ruta_raiz}pros/'
        self.ruta_runt_cpl = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/cpl/'
        self.ruta_data_wrst = ruta_wrst
        self.ruta_data_gfs = ruta_gfs
        self.name_gfs = gfs_file
        self.name_croco_his = f'{ruta_croco_his}/croco_his_{fecha_ini.strftime("%m%d")}.nc'
        self.name_ocenc = f'oce_{fecha_ini.strftime("%Y%m%d%H")}.nc'
        self.name_atmnc = f'atm_{fecha_ini.strftime("%Y%m%d%H")}.nc'
        self.base_croco = base_croco
        self.fini = fecha_ini
        self.name_atmsh = 'create_oasis_restart_from_preexisting_output_files.sh'
    def leer_his(self):
        nc_his = netCDF4.Dataset(self.name_croco_his,'r')
        time_his = nc_his.variables['time'][:]
        temp_his = nc_his.variables['temp'][:,:,:,:]
        mask_his = nc_his.variables['mask_rho'][:,:]
        lon_his = nc_his.variables['lon_rho'][:,:] + 360.
        lat_his = nc_his.variables['lat_rho'][:,:]
        date_croco = list(map(lambda x: self.base_croco + dt.timedelta(seconds=float(x)),time_his))
        nc_his.close()
        return temp_his, mask_his, lon_his, lat_his, date_croco
    def crear_oce(self):
        temp_his, mask_his, lon_his, lat_his, date_croco = self.leer_his()
        aux_dd = np.array(list(map(lambda x: x.day,date_croco)))
        aux_hh = np.array(list(map(lambda x: x.hour,date_croco)))
        pos_day = np.where((aux_dd==int(self.fini.day)) & (aux_hh==self.fini.hour))[0][0]
        sst_oce = temp_his[pos_day,-1,1:-1,1:-1]
        nw= netCDF4.Dataset(f'{self.ruta_runt_cpl}{self.name_ocenc}',mode='w',format='NETCDF4')
        dlon = nw.createDimension('eta_rho',sst_oce.shape[0])
        dlat = nw.createDimension('xi_rho', sst_oce.shape[1])
        data = nw.createVariable('SRMSSTV0','f4',('eta_rho','xi_rho'))
        data.long_name = 'potential temperature'
        data.units = 'Celsius'
        data.field = 'temperature, scalar, series'
        data.standard_name = 'sea_water_potential_temperature'
        data.coordinates = 'lat_rho lon_rho'
        data.cell_methods = 's_rho: mean time: mean'
        data[:,:] = sst_oce + 273.15
        nw.close()
        print ('Pos DATE = ',pos_day,date_croco[pos_day])
        print ('OCE = ',f'{self.ruta_runt_cpl}{self.name_ocenc}')
    def grafica_sst(self):
        nc_oce = netCDF4.Dataset(f'{self.ruta_runt_cpl}{self.name_ocenc}','r')
        sst = nc_oce.variables['SRMSSTV0'][:,:]
        nc_oce.close()
        temp_his, mask_his, lon_his, lat_his, date_croco = self.leer_his()

        fig = plt.figure(figsize=(10,6))
        ax = fig.add_subplot(111)
        mapa = ax.contourf(lon_his[1:-1,1:-1]-360.,lat_his[1:-1,1:-1],(sst-273.15)*mask_his[1:-1,1:-1],
                           np.linspace(22,34,26,endpoint='true'),cmap='jet')
        csmp = ax.contour(lon_his[1:-1,1:-1]-360.,lat_his[1:-1,1:-1],(sst-273.15)*mask_his[1:-1,1:-1],[25,27,29,31],colors='k')
        ax.clabel(csmp,inline=True,fmt='%2.1f',fontsize=10)
        plt.text(261.3-360,33.2,f'SST-oce.nc {self.fini.strftime("%Y %b %d %H:00")}',fontsize=20,fontweight='bold')
        cbar_ax = fig.add_axes([0.92, 0.14, 0.03, 0.68])
        cbar1 = plt.colorbar(mapa,cax=cbar_ax,orientation='vertical',ticks=[22,24,26,28,30,32])
        cbar1.set_label(r'SST $[C]$',fontsize=18)
        cbar1.ax.tick_params(labelsize=14)
        plt.savefig(f'{self.ruta_runt_cpl}oce_sst.png',dpi=200,bbox_inches="tight")
    def crear_atm(self):
        name_model = 'wrf'
        date_wrfrs = self.fini.strftime("%Y-%m-%d_%H:00:00")
        name_wrfrs = f'wrfrst_d01_{date_wrfrs}'
        com0 = f'{self.ruta_pros}{self.name_atmsh} {self.ruta_data_wrst}{name_wrfrs} {self.ruta_runt_cpl}{self.name_atmnc} {name_model}'
        os.system(com0)
        print (com0)

        grbs = pygrib.open(f'{self.ruta_data_gfs}{self.name_gfs}')
        aux_utau = grbs.select(name='Momentum flux, u component')[0]
        utau = aux_utau.values
        aux_vtau = grbs.select(name='Momentum flux, v component')[0]
        vtau = aux_vtau.values
        lat_mgfs, lon_mgfs = aux_vtau.latlons()
        grbs.close()
        print (f'*************** {name_wrfrs} ***************')
        nc_wrfi = netCDF4.Dataset(f'{self.ruta_data_wrst}{name_wrfrs}','r')
        latu = nc_wrfi.variables['XLAT_U'][0,:,:]
        lonu = nc_wrfi.variables['XLONG_U'][0,:,:]
        latv = nc_wrfi.variables['XLAT_V'][0,:,:]
        lonv = nc_wrfi.variables['XLONG_V'][0,:,:]
        xx_stag, yy_stag = np.meshgrid(lonu[0,:],latv[:,0])
        nc_wrfi.close()
        lon_var = lon_mgfs[0,:] - 360.
        lat_var = lat_mgfs[:,0]
        pos_lon = np.where((lon_var>np.min(lonu[0,:])) & (lon_var<np.max(lonu[0,:])))
        pos_lat = np.where((lat_var>np.min(latv[:,0])) & (lat_var<np.max(latv[:,0])))
        posx_cut, posy_cut = np.meshgrid(pos_lon[0],pos_lat[0])
        utau_cut = utau[posy_cut, posx_cut]
        vtau_cut = vtau[posy_cut, posx_cut]
        lon_cut = lon_var[pos_lon]
        lat_cut = lat_var[pos_lat]
        xx_cut, yy_cut = np.meshgrid(lon_cut,lat_cut)
        values_utau = np.reshape(utau_cut,-1)
        values_vtau = np.reshape(vtau_cut,-1)
        xx_vec = np.reshape(xx_cut,-1)
        yy_vec = np.reshape(yy_cut,-1)
        points = np.zeros((np.shape(xx_cut)[0]*np.shape(xx_cut)[1],2))
        points[:,0] = xx_vec
        points[:,1] = yy_vec
        utau_int = griddata(points, values_utau, (xx_stag, yy_stag), method='nearest')
        vtau_int = griddata(points, values_vtau, (xx_stag, yy_stag), method='nearest')

        nw= netCDF4.Dataset(f'{self.ruta_runt_cpl}{self.name_atmnc}',mode='a',format='NETCDF4')
        data = nw.createVariable('WRF_d01_EXT_d01_TAUX','f4',('south_north','west_east'))
        data.FieldType = 104.
        data.MemoryOrder = 'XY'
        data.description = 'ZONAL WIND STRESS'
        data.units = 'N m-2'
        data.stagger = ''
        data.coordinates = 'XLONG XLAT XTIME'
        data.cell_methods = 'Time: mean'
        data[:,:] = utau_int
        data = nw.createVariable('WRF_d01_EXT_d01_TAUY','f4',('south_north','west_east'))
        data.FieldType = 104.
        data.MemoryOrder = 'XY'
        data.description = 'MERIDIONAL WIND STRESS'
        data.units = 'N m-2'
        data.stagger = ''
        data.coordinates = 'XLONG XLAT XTIME'
        data.cell_methods = 'Time: mean'
        data[:,:] = vtau_int
        data = nw.createVariable('WRF_d01_EXT_d01_TAUMOD','f4',('south_north','west_east'))
        data.FieldType = 104.
        data.MemoryOrder = 'XY'
        data.description = 'WIND STRESS MODULE'
        data.units = 'N m-2'
        data.stagger = ''
        data.coordinates = 'XLONG XLAT XTIME'
        data.cell_methods = 'Time: mean'
        data[:,:] = np.zeros((356, 506))
        nw.close()

    def grafica_tau(self):
        nc_atm = netCDF4.Dataset(f'{self.ruta_runt_cpl}{self.name_atmnc}','r')
        taux_code = nc_atm.variables['WRF_d01_EXT_d01_TAUX'][:,:]
        tauy_code = nc_atm.variables['WRF_d01_EXT_d01_TAUY'][:,:]
        nc_atm.close()
        date_wrfrs = self.fini.strftime("%Y-%m-%d_%H:00:00")
        name_wrfrs = f'wrfrst_d01_{date_wrfrs}'
        nc_wrfi = netCDF4.Dataset(f'{self.ruta_data_wrst}{name_wrfrs}','r')
        latu = nc_wrfi.variables['XLAT_U'][0,:,:]
        lonu = nc_wrfi.variables['XLONG_U'][0,:,:]
        latv = nc_wrfi.variables['XLAT_V'][0,:,:]
        lonv = nc_wrfi.variables['XLONG_V'][0,:,:]
        xx_stag, yy_stag = np.meshgrid(lonu[0,:],latv[:,0])
        
        fig = plt.figure(figsize=(15,6))
        ax1 = fig.add_subplot(121)
        csx = ax1.contourf(xx_stag,yy_stag,taux_code,np.arange(-.5,.5,.001),cmap='PiYG',extend='both')
        ax1.set_title('TauX',fontsize=20)
        ax2 = fig.add_subplot(122)
        csy = ax2.contourf(xx_stag,yy_stag,tauy_code,np.arange(-.5,.5,.001),cmap='bwr',extend='both')
        ax2.set_title('TauY',fontsize=20)
        cbar_ax = fig.add_axes([.92, 0.16, 0.025, 0.62])
        cbar = plt.colorbar(csx,cax=cbar_ax,orientation='vertical')
        cbar.set_label(r'Taux $[Nm^{-2}]$',fontsize=18)
        cbar_ax = fig.add_axes([1.01, 0.16, 0.025, 0.62])
        cbar = plt.colorbar(csy,cax=cbar_ax,orientation='vertical')
        cbar.set_label(r'Tauy $[Nm^{-2}]$',fontsize=18)
        plt.savefig(f'{self.ruta_runt_cpl}atm_tau.png',dpi=200,bbox_inches="tight")

class pre_croco():
    def __init__ (self,ruta_raiz,fecha_ini,ruta_croco_his,base_croco):
        self.ruta_runt_cpl = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/cpl/'
        self.name_croco_his = f'{ruta_croco_his}/croco_his_{fecha_ini.strftime("%m%d")}.nc'
        self.name_croco_ini = f'croco_ini_{fecha_ini.strftime("%Y%m%d%H")}.nc'
        self.name_croco_bry = f'croco_bry_{fecha_ini.strftime("%Y%m%d%H")}.nc'
        self.base_croco = base_croco
        self.fini = fecha_ini
    def leer_his(self):
        nc_his = netCDF4.Dataset(self.name_croco_his,'r')
        time_his = nc_his.variables['time'][:]
        temp_his = nc_his.variables['temp'][:,:,:,:]
        mask_his = nc_his.variables['mask_rho'][:,:]
        lon_his = nc_his.variables['lon_rho'][:,:] + 360.
        lat_his = nc_his.variables['lat_rho'][:,:]
        date_croco = list(map(lambda x: self.base_croco + dt.timedelta(seconds=float(x)),time_his))
        nc_his.close()
        return temp_his, mask_his, lon_his, lat_his, date_croco
    def crear_ini(self):
        nc = netCDF4.Dataset(self.name_croco_his,'r')
        temp_his = nc.variables['temp'][:,:,:,:]
        salt_his = nc.variables['salt'][:,:,:,:]
        uu_his = nc.variables['u'][:,:,:,:]
        vv_his = nc.variables['v'][:,:,:,:]
        ub_his = nc.variables['ubar'][:,:,:]
        vb_his = nc.variables['vbar'][:,:,:]
        zeta_his = nc.variables['zeta'][:,:,:]
        time_his = nc.variables['time'][:]
        csr_his = nc.variables['Cs_r'][:]
        scr_his = nc.variables['sc_r'][:]
        nc.close()

        date_croco = list(map(lambda x: self.base_croco + dt.timedelta(seconds=float(x)),time_his))
        aux_dd = np.array(list(map(lambda x: x.day,date_croco)))
        aux_hh = np.array(list(map(lambda x: x.hour,date_croco)))
        pos_day = np.where((aux_dd==int(self.fini.day)) & (aux_hh==self.fini.hour))[0][0]

        nw = netCDF4.Dataset(f'{self.ruta_runt_cpl}croco_files/{self.name_croco_ini}',mode='w',format='NETCDF4')
        d1 = nw.createDimension('xi_u',241)
        d2 = nw.createDimension('xi_v',242)
        d3 = nw.createDimension('xi_rho',242)
        d4 = nw.createDimension('eta_u',173)
        d5 = nw.createDimension('eta_v',172)
        d6 = nw.createDimension('eta_rho',173)
        d7 = nw.createDimension('s_rho',32)
        d8 = nw.createDimension('s_w',33)
        d9 = nw.createDimension('tracer',2)
        dm = nw.createDimension('time',1)
        d0 = nw.createDimension('one',1)
        data = nw.createVariable('spherical','S1',('one'))
        data ='T'
        data = nw.createVariable('Vtransform','i4',('one'))
        data.long_name = 'vertical terrain-following transformation equation'
        data = 2
        data = nw.createVariable('Vstretching','i4',('one'))
        data.long_name = 'vertical terrain-following stretching function'
        data = 1
        data = nw.createVariable('tstart','f8',('one'))
        data.long_name = 'start processing day'
        data.units = 'day'
        data[:] = 0.
        data = nw.createVariable('tend','f8',('one'))
        data.long_name = 'end processing day'
        data.units = 'day'
        data[:] = 0.
        data = nw.createVariable('theta_s','f8',('one'))
        data.long_name = 'S-coordinate surface control parameter'
        data.units = 'nondimensional'
        data[:] = 7.
        data = nw.createVariable('theta_b','f8',('one'))
        data.long_name = 'S-coordinate bottom control parameter'
        data.units = 'nondimensional'
        data[:] = 2.
        data = nw.createVariable('Tcline','f8',('one'))
        data.long_name = 'S-coordinate surface/bottom layer width'
        data.units = 'meter'
        data[:] = 200.
        data = nw.createVariable('hc','f8',('one'))
        data.long_name = 'S-coordinate parameter, critical depth'
        data.units = 'meter'
        data[:] = 200.
        data = nw.createVariable('sc_r','f8',('s_rho'))
        data.long_name = 'S-coordinate at RHO-points'
        data.units = 'nondimensional'
        data.valid_min = -1.
        data.valid_max = 0.
        data[:] = scr_his
        data = nw.createVariable('Cs_r','f8',('s_rho'))
        data.long_name = 'S-coordinate stretching curves at RHO-points'
        data.units = 'nondimensional'
        data.valid_min = -1.
        data.valid_max = 0.
        data[:] = csr_his
        data = nw.createVariable('ocean_time','f8',('time'))
        data.long_name = 'time since initialization'
        data.units = 'second'
        data[:] = 0.
        data = nw.createVariable('scrum_time','f8',('time'))
        data.long_name = 'time since initialization'
        data.units = 'second'
        data[:] = 0.
        data = nw.createVariable('u','f8',('time','s_rho','eta_u','xi_u'))
        data.long_name = 'u-momentum component'
        data.units = 'meter second-1'
        data[0,:,:,:] = uu_his[pos_day,:,:,:]
        data = nw.createVariable('v','f8',('time','s_rho','eta_v','xi_v'))
        data.long_name = 'v-momentum component'
        data.units = 'meter second-1'
        data[0,:,:,:] = vv_his[pos_day,:,:,:]
        data = nw.createVariable('ubar','f8',('time','eta_u','xi_u'))
        data.long_name = 'vertically integrated u-momentum component'
        data.units = 'meter second-1'
        data[0,:,:] = ub_his[pos_day,:,:]
        data = nw.createVariable('vbar','f8',('time','eta_v','xi_v'))
        data.long_name = 'vertically integrated v-momentum component'
        data.units = 'meter second-1'
        data[0,:,:] = vb_his[pos_day,:,:]
        data = nw.createVariable('zeta','f8',('time','eta_rho','xi_rho'))
        data.long_name = 'free-surface'
        data.units = 'meter'
        data[0,:,:] = zeta_his[pos_day,:,:]
        data = nw.createVariable('temp','f8',('time','s_rho','eta_rho','xi_rho'))
        data.long_name = 'potential temperature'
        data.units = 'Celsius'
        data[:,:,:,:] = temp_his[pos_day,:,:,:]
        data = nw.createVariable('salt','f8',('time','s_rho','eta_rho','xi_rho'))
        data.long_name = 'salinity'
        data.units = 'PSU'
        data[0,:,:,:] = salt_his[pos_day,:,:,:]
        nw.close()
        print ('INI = ',pos_day,date_croco[pos_day])
        print ('INI = ',f'{self.ruta_runt_cpl}{self.name_croco_ini}')

    def grafica_ini(self):
        nc_ini = netCDF4.Dataset(f'{self.ruta_runt_cpl}croco_files/{self.name_croco_ini}','r')
        temp_ini = nc_ini.variables['temp'][:,:,:,:]
        nc_ini.close()
        temp_his, mask_his, lon_his, lat_his, date_croco = self.leer_his()
        fig = plt.figure(figsize=(10,6))
        ax = fig.add_subplot(111)
        mapa = ax.contourf(lon_his-360.,lat_his,(temp_ini[0,-1,:,:])*mask_his,
                           np.linspace(22,34,26,endpoint='true'),cmap='jet')
        csmp = ax.contour(lon_his-360.,lat_his,(temp_ini[0,-1,:,:])*mask_his,[25,27,29,31],colors='k')
        ax.clabel(csmp,inline=True,fmt='%2.1f',fontsize=10)
        plt.text(261.3-360,33.2,f'SST-croco_ini.nc {self.fini.strftime("%Y %b %d %H:00")}',fontsize=20,fontweight='bold')
        cbar_ax = fig.add_axes([0.92, 0.14, 0.03, 0.68])
        cbar1 = plt.colorbar(mapa,cax=cbar_ax,orientation='vertical',ticks=[22,24,26,28,30,32])
        cbar1.set_label(r'SST $[C]$',fontsize=18)
        cbar1.ax.tick_params(labelsize=14)
        plt.savefig(f'{self.ruta_runt_cpl}croco_ini.png',dpi=200,bbox_inches="tight")

    def crear_bry(self):
        time_aux = np.array((0,1),dtype='float')
        date_bry = np.array(list(map(lambda x: self.fini + dt.timedelta(days=x),time_aux)))
        temp_aux = np.zeros((len(time_aux),32,173,242))
        salt_aux = np.zeros((len(time_aux),32,173,242))
        uu_aux = np.zeros((len(time_aux),32,173,241))
        vv_aux = np.zeros((len(time_aux),32,172,242))
        ub_aux = np.zeros((len(time_aux),173,241))
        vb_aux = np.zeros((len(time_aux),172,242))
        zeta_aux = np.zeros((len(time_aux),173,242))
        csr_aux = np.zeros((len(time_aux),32))
        scr_aux = np.zeros((len(time_aux),32))
        scw_aux = np.zeros((len(time_aux),33))
        csw_aux = np.zeros((len(time_aux),33))

        nc = netCDF4.Dataset(self.name_croco_his,'r')
        time_his = nc.variables['time'][:]

        date_croco = np.array(list(map(lambda x: self.base_croco + dt.timedelta(seconds=float(x)),time_his)))
        aux_dd = np.array(list(map(lambda x: x.day,date_croco)))
        aux_hh = np.array(list(map(lambda x: x.hour,date_croco)))

        for jj in range(len(date_bry)):
            pos_day = np.where((aux_dd==int(date_bry[jj].day)) & (aux_hh==date_bry[jj].hour))[0][0]
            print ('BRY ',date_croco[pos_day])
            temp_aux[jj,:,:,:] = nc.variables['temp'][pos_day,:,:,:]
            salt_aux[jj,:,:,:] = nc.variables['salt'][pos_day,:,:,:]
            uu_aux[jj,:,:,:] = nc.variables['u'][pos_day,:,:,:]
            vv_aux[jj,:,:,:] = nc.variables['v'][pos_day,:,:,:]
            ub_aux[jj,:,:] = nc.variables['ubar'][pos_day,:,:]
            vb_aux[jj,:,:] = nc.variables['vbar'][pos_day,:,:]
            zeta_aux[jj,:,:] = nc.variables['zeta'][pos_day,:,:]
            csr_aux[jj,:] = nc.variables['Cs_r'][:]
            scr_aux[jj,:] = nc.variables['sc_r'][:]
            scw_aux[jj,:] = nc.variables['sc_w'][:]
            csw_aux[jj,:] = nc.variables['Cs_w'][:]
        nc.close()

        nw = netCDF4.Dataset(f'{self.ruta_runt_cpl}croco_files/{self.name_croco_bry}',mode='w',format='NETCDF4')
        d1 = nw.createDimension('xi_u',241)
        d2 = nw.createDimension('xi_v',242)
        d3 = nw.createDimension('xi_rho',242)
        d4 = nw.createDimension('eta_u',173)
        d5 = nw.createDimension('eta_v',172)
        d6 = nw.createDimension('eta_rho',173)
        d7 = nw.createDimension('s_rho',32)
        d8 = nw.createDimension('s_w',33)
        d9 = nw.createDimension('tracer',2)
        d0 = nw.createDimension('one',1)
        tbry = nw.createDimension('bry_time',len(time_aux))
        tclm = nw.createDimension('tclm_time',len(time_aux)) 
        ttem = nw.createDimension('temp_time',len(time_aux))
        tslm = nw.createDimension('sclm_time',len(time_aux))
        tsal = nw.createDimension('salt_time',len(time_aux)) 
        tucl = nw.createDimension('uclm_time',len(time_aux)) 
        tvcl = nw.createDimension('vclm_time',len(time_aux)) 
        tv2d = nw.createDimension('v2d_time',len(time_aux)) 
        tv3d = nw.createDimension('v3d_time',len(time_aux)) 
        tssh = nw.createDimension('ssh_time',len(time_aux)) 
        tzet = nw.createDimension('zeta_time',len(time_aux)) 

        data = nw.createVariable('spherical','S1',('one'))
        data.long_name = 'grid type logical switch'
        data.flag_values = 'T, F'
        data.flag_meanings = 'spherical Cartesian'
        data ='T'
        data = nw.createVariable('Vtransform','i4',('one'))
        data.long_name = 'vertical terrain-following transformation equation'
        data = 2
        data = nw.createVariable('Vstretching','i4',('one'))
        data.long_name = 'vertical terrain-following stretching function'
        data = 1
        data = nw.createVariable('tstart','f8',('one'))
        data.long_name = 'start processing day'
        data.units = 'day'
        data[:] = 0.
        data = nw.createVariable('tend','f8',('one'))
        data.long_name = 'end processing day'
        data.units = 'day'
        data[:] = 10.
        data = nw.createVariable('theta_s','f8',('one'))
        data.long_name = 'S-coordinate surface control parameter'
        data.units = 'nondimensional'
        data[:] = 7.
        data = nw.createVariable('theta_b','f8',('one'))
        data.long_name = 'S-coordinate bottom control parameter'
        data.units = 'nondimensional'
        data[:] = 2.
        data = nw.createVariable('Tcline','f8',('one'))
        data.long_name = 'S-coordinate surface/bottom layer width'
        data.units = 'meter'
        data[:] = 200.
        data = nw.createVariable('hc','f8',('one'))
        data.long_name = 'S-coordinate parameter, critical depth'
        data.units = 'meter'
        data[:] = 200.
        data = nw.createVariable('sc_r','f8',('s_rho'))
        data.long_name = 'S-coordinate at RHO-points'
        data.units = 'nondimensional'
        data.valid_min = -1.
        data.valid_max = 0.
        data.positive = 'up'
        data.standard_name = 'ocena_s_coordinate_g2'
        data.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
        data[:] = scr_aux[0,:] 
        data = nw.createVariable('sc_w','f8',('s_w'))
        data.long_name = 'S-coordinate at W-points'
        data.units = 'nondimensional'
        data.valid_min = -1.
        data.valid_max = 0.
        data.positive = 'up'
        data.standard_name = 'ocena_s_coordinate_g2'
        data.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
        data[:] = scw_aux[0,:]
        data = nw.createVariable('Cs_r','f8',('s_rho'))
        data.long_name = 'S-coordinate stretching curves at RHO-points'
        data.units = 'nondimensional'
        data.valid_min = -1.
        data.valid_max = 0.
        data[:] = csr_aux[0,:]
        data = nw.createVariable('Cs_w','f8',('s_w'))
        data.long_name = 'S-coordinate stretching curves at W-points'
        data.units = 'nondimensional'
        data.valid_min = -1.
        data.valid_max = 0.
        data[:] = csw_aux[0,:]
        data = nw.createVariable('bry_time','f8',('bry_time'))
        data.long_name = 'time for boundary climatology'
        data.units = 'day'
        data.calendar = '360.0 days in every year'
        data.cycle_length = 10.
        data[:] = time_aux
        data = nw.createVariable('tclm_time','f8',('tclm_time'))
        data.long_name = 'time for temperature climatology'
        data.units = 'day'
        data.calendar = '360.0 days in every year'
        data.cycle_length = 10.
        data[:] = time_aux
        data = nw.createVariable('temp_time','f8',('temp_time'))
        data.long_name = 'time for temperature climatology'
        data.units = 'day'
        data.calendar = '360.0 days in every year'
        data.cycle_length = 10.
        data[:] = time_aux
        data = nw.createVariable('sclm_time','f8',('sclm_time'))
        data.long_name = 'time for salinity climatology'
        data.units = 'day'
        data.calendar = '360.0 days in every year'
        data.cycle_length = 10.
        data[:] = time_aux
        data = nw.createVariable('salt_time','f8',('salt_time'))
        data.long_name = 'time for salinity climatology'
        data.units = 'day'
        data.calendar = '360.0 days in every year'
        data.cycle_length = 10.
        data[:] = time_aux
        data = nw.createVariable('uclm_time','f8',('uclm_time'))
        data.long_name = 'time climatological u'
        data.units = 'day'
        data.calendar = '360.0 days in every year'
        data.cycle_length = 10.
        data[:] = time_aux
        data = nw.createVariable('vclm_time','f8',('vclm_time'))
        data.long_name = 'time climatological v'
        data.units = 'day'
        data.calendar = '360.0 days in every year'
        data.cycle_length = 10.
        data[:] = time_aux
        data = nw.createVariable('v2d_time','f8',('v2d_time'))
        data.long_name = 'time for 2D velocity climatology'
        data.units = 'day'
        data.calendar = '360.0 days in every year'
        data.cycle_length = 10.
        data[:] = time_aux
        data = nw.createVariable('v3d_time','f8',('v3d_time'))
        data.long_name = 'time for 3D velocity climatology'
        data.units = 'day'
        data.calendar = '360.0 days in every year'
        data.cycle_length = 10.
        data[:] = time_aux
        data = nw.createVariable('ssh_time','f8',('ssh_time'))
        data.long_name = 'time for sea surface height'
        data.units = 'day'
        data.calendar = '360.0 days in every year'
        data.cycle_length = 10.
        data[:] = time_aux
        data = nw.createVariable('zeta_time','f8',('zeta_time'))
        data.long_name = 'time for sea surface height'
        data.units = 'day'
        data.calendar = '360.0 days in every year'
        data.cycle_length = 10.
        data[:] = time_aux
        data = nw.createVariable('temp_east','f8',('temp_time','s_rho','eta_rho'))
        data.long_name = 'eastern boundary potential temperature' 
        data.units = 'Celsius'
        data.coordinates = 'lat_rho s_rho temp_time'
        for jj in range(len(time_aux)):
            data[jj,:,:] = temp_aux[jj,:,:,-1]

        data = nw.createVariable('salt_east','f8',('salt_time','s_rho','eta_rho'))
        data.long_name = 'eastern boundary salinity' 
        data.units = 'PSU'
        data.coordinates = 'lat_rho s_rho salt_time'
        for jj in range(len(time_aux)):
            data[jj,:,:] = salt_aux[jj,:,:,-1]

        data = nw.createVariable('u_east','f8',('v3d_time','s_rho','eta_rho'))
        data.long_name = 'eastern boundary u-momentum component' 
        data.units = 'meter second-1'
        data.coordinates = 'lat_u s_rho u_time'
        for jj in range(len(time_aux)):
            data[jj,:,:] = uu_aux[jj,:,:,-1]

        data = nw.createVariable('v_east','f8',('v3d_time','s_rho','eta_v'))
        data.long_name = 'eastern boundary v-momentum component' 
        data.units = 'meter second-1'
        data.coordinates = 'lat_v s_rho vclm_time'
        for jj in range(len(time_aux)):
            data[jj,:,:] = vv_aux[jj,:,:,-1]

        data = nw.createVariable('ubar_east','f8',('v2d_time','eta_rho'))
        data.long_name = 'eastern boundary vertically integrated u-momentum component' 
        data.units = 'meter second-1'
        data.coordinates = 'lat_u uclm_time'
        for jj in range(len(time_aux)):
            data[jj,:] = ub_aux[jj,:,-1]

        data = nw.createVariable('vbar_east','f8',('v2d_time','eta_v'))
        data.long_name = 'eastern boundary vertically integrated v-momentum component' 
        data.units = 'meter second-1'
        data.coordinates = 'lat_v vclm_time'
        for jj in range(len(time_aux)):
            data[jj,:] = vb_aux[jj,:,-1]

        data = nw.createVariable('zeta_east','f8',('zeta_time','eta_rho'))
        data.long_name = 'eastern boundary sea surface height' 
        data.units = 'meter'
        data.coordinates = 'lat_rho zeta_time'
        for jj in range(len(time_aux)):
            data[jj,:] = zeta_aux[jj,:,-1]

        data = nw.createVariable('temp_north','f8',('temp_time','s_rho','xi_rho'))
        data.long_name = 'northern boundary potential temperature' 
        data.units = 'Celsius'
        data.coordinates = 'lon_rho s_rho temp_time'
        for jj in range(len(time_aux)):
            data[jj,:,:] = temp_aux[jj,:,-1,:]

        data = nw.createVariable('salt_north','f8',('salt_time','s_rho','xi_rho'))
        data.long_name = 'northern boundary salinity' 
        data.units = 'PSU'
        data.coordinates = 'lon_rho s_rho salt_time'
        for jj in range(len(time_aux)):
            data[jj,:,:] = salt_aux[jj,:,-1,:]

        data = nw.createVariable('u_north','f8',('v3d_time','s_rho','xi_u'))
        data.long_name = 'northern boundary u-momentum component' 
        data.units = 'meter second-1'
        data.coordinates = 'lon_u s_rho u_time'
        for jj in range(len(time_aux)):
            data[jj,:,:] = uu_aux[jj,:,-1,:]

        data = nw.createVariable('v_north','f8',('v3d_time','s_rho','xi_rho'))
        data.long_name = 'northern boundary v-momentum component' 
        data.units = 'meter second-1'
        data.coordinates = 'lon_v s_rho vclm_time'
        for jj in range(len(time_aux)):
            data[jj,:,:] = vv_aux[jj,:,-1,:]

        data = nw.createVariable('ubar_north','f8',('v2d_time','xi_u'))
        data.long_name = 'northern boundary vertically integrated u-momentum component' 
        data.units = 'meter second-1'
        data.coordinates = 'lon_u uclm_time'
        for jj in range(len(time_aux)):
            data[jj,:] = ub_aux[jj,-1,:]

        data = nw.createVariable('vbar_north','f8',('v2d_time','xi_rho'))
        data.long_name = 'northern boundary vertically integrated v-momentum component' 
        data.units = 'meter second-1'
        data.coordinates = 'lat_v vclm_time'
        for jj in range(len(time_aux)):
            data[jj,:] = vb_aux[jj,-1,:]

        data = nw.createVariable('zeta_north','f8',('zeta_time','xi_rho'))
        data.long_name = 'northern boundary sea surface height' 
        data.units = 'meter'
        data.coordinates = 'lon_rho zeta_time'
        for jj in range(len(time_aux)):
            data[jj,:] = zeta_aux[jj,-1,:]
        nw.close()

        print ('INI = ',f'{self.ruta_runt_cpl}{self.name_croco_bry}')

    def grafica_bry(self):
        time_aux = np.array((0,1),dtype='float')
        nc = netCDF4.Dataset(f'{self.ruta_runt_cpl}croco_files/{self.name_croco_bry}',mode='r',format='NETCDF4')
        temp_east = nc.variables['temp_east'][:,:,:]
        temp_nort = nc.variables['temp_north'][:,:,:]
        salt_east = nc.variables['salt_east'][:,:,:]
        salt_nort = nc.variables['salt_north'][:,:,:]
        v_nort = nc.variables['v_north'][:,:,:]
        v_east = nc.variables['v_east'][:,:,:]
        u_nort = nc.variables['u_north'][:,:,:]
        u_east = nc.variables['u_east'][:,:,:]
        nc.close()

        fig = plt.figure(figsize=(10,6))
        for jj in range(len(time_aux)):
            ax = fig.add_subplot(1,2,jj+1)
            cs = ax.contourf(temp_nort[jj,:,:],np.arange(0,32,1),cmap='jet')
            ax.contour(temp_nort[jj,:,:],[20],colors='k')
        cbar_ax = fig.add_axes([0.92, 0.125, 0.04, 0.755])
        cbar = plt.colorbar(cs,cax=cbar_ax,orientation='vertical')
        cbar.set_label(r'Temperature $[^{\circ}C]$',fontsize=22)
        plt.savefig(f'{self.ruta_runt_cpl}croco_bry_temp.png',dpi=200,bbox_inches="tight")

        fig = plt.figure(figsize=(10,6))
        for jj in range(len(time_aux)):
            ax = fig.add_subplot(1,2,jj+1)
            cs = ax.contourf(salt_nort[jj,:,:],np.arange(30,37.5,.1),cmap='brg')
            ax.contour(salt_nort[jj,:,:],[20],colors='k')
        cbar_ax = fig.add_axes([0.92, 0.125, 0.04, 0.755])
        cbar = plt.colorbar(cs,cax=cbar_ax,orientation='vertical')
        cbar.set_label(r'Salinity $[psu]$',fontsize=22)
        plt.savefig(f'{self.ruta_runt_cpl}croco_bry_salt.png',dpi=200,bbox_inches="tight")


class ncversion():
    def __init__ (self,ruta_raiz,fecha_ini):
        self.ruta_runt_cpl = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/cpl/'
        self.name_croco_ini = f'croco_ini_{fecha_ini.strftime("%Y%m%d%H")}.nc'
        self.name_croco_bry = f'croco_bry_{fecha_ini.strftime("%Y%m%d%H")}.nc'
        self.name_croco_grd = 'croco_grd.nc'
        self.name_oce = f'oce_{fecha_ini.strftime("%Y%m%d%H")}.nc'
        self.name_atm = f'atm_{fecha_ini.strftime("%Y%m%d%H")}.nc'
        self.wrfinput = 'wrfinput_d01'
        self.wrflowinp = 'wrflowinp_d01'
        self.wrfbdy = 'wrfbdy_d01'
        self.fini = fecha_ini
    def mover_nc(self):
        subprocess.call(['mkdir',f'{self.ruta_runt_cpl}nc_orig/'])
        subprocess.call(['mv',f'{self.ruta_runt_cpl}croco_files/{self.name_croco_ini}',f'{self.ruta_runt_cpl}nc_orig/'])
        subprocess.call(['mv',f'{self.ruta_runt_cpl}croco_files/{self.name_croco_bry}',f'{self.ruta_runt_cpl}nc_orig/'])
        subprocess.call(['mv',f'{self.ruta_runt_cpl}croco_files/{self.name_croco_grd}',f'{self.ruta_runt_cpl}nc_orig/'])
        subprocess.call(['mv',f'{self.ruta_runt_cpl}{self.name_oce}',f'{self.ruta_runt_cpl}nc_orig/'])
        subprocess.call(['mv',f'{self.ruta_runt_cpl}{self.name_atm}',f'{self.ruta_runt_cpl}nc_orig/'])
        subprocess.call(['mv',f'{self.ruta_runt_cpl}{self.wrfinput}',f'{self.ruta_runt_cpl}nc_orig/'])
        subprocess.call(['mv',f'{self.ruta_runt_cpl}{self.wrflowinp}',f'{self.ruta_runt_cpl}nc_orig/'])
        subprocess.call(['mv',f'{self.ruta_runt_cpl}{self.wrfbdy}',f'{self.ruta_runt_cpl}nc_orig/'])
    def netcdf3(self):
        archivos = os.listdir(f'{self.ruta_runt_cpl}nc_orig/')
        for archivo in archivos:
            ruta_archivo = os.path.join(f'{self.ruta_runt_cpl}nc_orig/', archivo)
            if ruta_archivo.split('/')[-1][0:3] == 'atm':
                name_out = f'{self.ruta_runt_cpl}atm.nc'
                com0 = f'nccopy -k 1 {ruta_archivo} {name_out}'
                os.system(com0)
                print (name_out)
            if ruta_archivo.split('/')[-1][0:3] == 'oce':
                name_out = f'{self.ruta_runt_cpl}oce.nc'
                com0 = f'nccopy -k 1 {ruta_archivo} {name_out}'
                os.system(com0)
                print (name_out)
            if ruta_archivo.split('/')[-1][6:9] == 'grd':
                name_out = f'{self.ruta_runt_cpl}croco_files/croco_grd.nc'
                com0 = f'nccopy -k 1 {ruta_archivo} {name_out}'
                os.system(com0)
                print (name_out)
            if ruta_archivo.split('/')[-1][6:9] == 'ini':
                name_out = f'{self.ruta_runt_cpl}croco_files/croco_ini.nc'
                com0 = f'nccopy -k 1 {ruta_archivo} {name_out}'
                os.system(com0)
                print (name_out)
            if ruta_archivo.split('/')[-1][6:9] == 'bry':
                name_out = f'{self.ruta_runt_cpl}croco_files/croco_bry.nc'
                com0 = f'nccopy -k 1 {ruta_archivo} {name_out}'
                os.system(com0)
                print (name_out)
            if ruta_archivo.split('/')[-1][0:9] == 'wrflowinp':
                name_out = f'{self.ruta_runt_cpl}wrflowinp_d01'
                com0 = f'nccopy -k 1 {ruta_archivo} {name_out}'
                os.system(com0)
                print (name_out)
            if ruta_archivo.split('/')[-1][0:8] == 'wrfinput':
                name_out = f'{self.ruta_runt_cpl}wrfinput_d01'
                com0 = f'nccopy -k 1 {ruta_archivo} {name_out}'
                os.system(com0)
                print (name_out)
            if ruta_archivo.split('/')[-1][0:6] == 'wrfbdy':
                name_out = f'{self.ruta_runt_cpl}wrfbdy_d01'
                com0 = f'nccopy -k 1 {ruta_archivo} {name_out}'
                os.system(com0)
                print (name_out)

class cpl_wrf():
    def __init__ (self,ruta_raiz,fecha_ini,nprocw,nprocc):
        self.ruta_runt_cpl = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/cpl/'
        self.nprocw = nprocw
        self.nprocc = nprocc
    def ejecutar_wrf(self):
        com0 = f'mpirun -np {self.nprocw} {self.ruta_runt_cpl}wrfexe  : -np {self.nprocc} {self.ruta_runt_cpl}crocox > {self.ruta_runt_cpl}coupled.log'
        os.chdir(self.ruta_runt_cpl)
        os.system(com0)

class oper_resultados():
    def __init__ (self,ruta_raiz,fecha_ini):
        self.ruta_runt_cpl = f'{ruta_raiz}run/{fecha_ini.strftime("%Y%m%d%H")}/cpl/'
        self.ruta_data_cpl = f'{ruta_raiz}data/{fecha_ini.strftime("%Y%m%d%H")}/cpl/'
    def copia_wrfout(self):
        archivos = os.listdir(f'{self.ruta_runt_cpl}')
        for archivo in archivos:
            ruta_archivo = os.path.join(self.ruta_runt_cpl, archivo)
            if os.path.isfile(ruta_archivo):
                if ruta_archivo.split('/')[-1][0:6] == 'wrfout':
                    shutil.copy(ruta_archivo, self.ruta_data_cpl)
    def copia_wrfhr(self):
        archivos = os.listdir(f'{self.ruta_runt_cpl}')
        for archivo in archivos:
            ruta_archivo = os.path.join(self.ruta_runt_cpl, archivo)
            if os.path.isfile(ruta_archivo):
                if ruta_archivo.split('/')[-1][0:5] == 'wrfhr':
                    shutil.copy(ruta_archivo, self.ruta_data_cpl)
    def copia_croco(self):
        archivos = os.listdir(f'{self.ruta_runt_cpl}croco_outfl/')
        print (archivos, f'{self.ruta_runt_cpl}croco_outfl/')
        for archivo in archivos:
            ruta_archivo = os.path.join(f'{self.ruta_runt_cpl}croco_outfl/', archivo)
            if os.path.isfile(ruta_archivo):
                print (ruta_archivo)
                shutil.copy(ruta_archivo, self.ruta_data_cpl)
    def copia_config(self):
        archivos = os.listdir(f'{self.ruta_runt_cpl}')
        for archivo in archivos:
            ruta_archivo = os.path.join(self.ruta_runt_cpl, archivo)
            if os.path.isfile(ruta_archivo):
                if(ruta_archivo.split('/')[-1] == 'namelist.input'):
                    shutil.copy(ruta_archivo, self.ruta_data_cpl)



