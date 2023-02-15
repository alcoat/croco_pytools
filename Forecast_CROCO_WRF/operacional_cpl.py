# -*- coding: utf-8 -*-
#!/usr/bin/python

import datetime as dt
#import clases_wrf as wrf_oper
import clases_cpl as cpl_oper


# ============================= Users changes ===================

ruta_main = '/exports/home/mzapatahe/piloto/'
ruta_croco_his = '/exports/home/mzapatahe/piloto/data/2022092112/croco/'

nhora = 72
fecha_ini = dt.datetime(2022, 9, 21, 12, 0)
fecha_fin = fecha_ini + dt.timedelta(hours=nhora)
fecha_rst = fecha_ini - dt.timedelta(hours=24)
base_croco = dt.datetime(2005,1,1,0)

# ===============================================================

copiar = cpl_oper.oper_copiar(ruta_raiz=ruta_main,
                              fecha_ini=fecha_ini)

copiar.copia_wrf()
copiar.copia_croco()
copiar.copia_oasis()

prewrf = cpl_oper.pre_wrf(ruta_raiz=ruta_main,
                          fecha_ini=fecha_ini,
                          fecha_fin=fecha_fin)
prewrf.copia_wrf_archivos()
prewrf.editar_namelist()
prewrf.modifica_wrfinput()
prewrf.grafica_mascara()

preoasis = cpl_oper.pre_oasis(ruta_raiz=ruta_main,
                              fecha_ini=fecha_ini,
                              ruta_croco_his=ruta_croco_his,
                              base_croco=base_croco,
                              ruta_wrst=f'{ruta_main}data/{fecha_rst.strftime("%Y%m%d%H")}/wrf/',
                              ruta_gfs=f'{ruta_main}data/{fecha_rst.strftime("%Y%m%d%H")}/gfs/',
                              gfs_file='gfs.t12z.pgrb2.0p25.f024')

preoasis.crear_oce()
preoasis.grafica_sst()
preoasis.crear_atm()
preoasis.grafica_tau()

precroco = cpl_oper.pre_croco(ruta_raiz=ruta_main,
                              fecha_ini=fecha_ini,
                              ruta_croco_his=ruta_croco_his,
                              base_croco=base_croco)
precroco.crear_ini()
precroco.grafica_ini()
precroco.crear_bry()
precroco.grafica_bry()


version = cpl_oper.ncversion(ruta_raiz=ruta_main,
                             fecha_ini=fecha_ini)
version.mover_nc()
version.netcdf3()

cpl = cpl_oper.cpl_wrf(ruta_raiz=ruta_main,
                       fecha_ini=fecha_ini,
                       nprocw='56',
                       nprocc='8')
cpl.ejecutar_wrf()

copiar = cpl_oper.oper_resultados(ruta_raiz=ruta_main,
                                  fecha_ini=fecha_ini)

copiar.copia_wrfout()
copiar.copia_wrfhr()
copiar.copia_croco()
copiar.copia_config()
