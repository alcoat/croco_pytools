# -*- coding: utf-8 -*-
#!/usr/bin/python

import datetime as dt
import clases_wrf as wrf_oper

# ============================= Users changes ===================

ruta_main = '/exports/home/ffayalac/runs/ANH/piloto/'
ruta_wpsx = '/opt/ohpc/WRF/4.1.5/WPS/'
ruta_arwx = '/opt/ohpc/WRF/4.2.1/WRF/main/'
ruta_mpir = '/usr/lib64/mpich-3.2/bin/mpirun'

nhora = 72
fecha_ini = dt.datetime(2022, 11, 22, 12, 0)
fecha_fin = fecha_ini + dt.timedelta(hours=nhora)

# ===============================================================

carpeteador = wrf_oper.oper_carpetas(ruta_raiz=ruta_main,
                                     fecha_ini=fecha_ini)

carpeteador.crear_carpetas_l1()
carpeteador.crear_carpetas_l2()
carpeteador.crear_carpetas_l3()

gfs = wrf_oper.wrf_gfs(ruta_raiz=ruta_main,
                       fecha_ini=fecha_ini,
                       nh=nhora)
gfs.descarga_gfs()


geogrid = wrf_oper.wrf_geo(ruta_raiz=ruta_main,
                           ruta_geox=ruta_wpsx,
                           fecha_ini=fecha_ini)
geogrid.copiar_archivos()
geogrid.editar_namelist()
geogrid.ejecutar_geogrid()

linkgrib = wrf_oper.wrf_lng(ruta_raiz=ruta_main,
                            ruta_lngx=ruta_wpsx,
                            fecha_ini=fecha_ini,
                            suffix='gfs')
linkgrib.ejecutar_linkgrib()

ungrib = wrf_oper.wrf_ung(ruta_raiz=ruta_main,
                          ruta_ungx=ruta_wpsx,
                          fecha_ini=fecha_ini,
                          fecha_fin=fecha_fin)
ungrib.copiar_archivos()
ungrib.editar_namelist()
ungrib.ejecutar_ungrib()

metgrid = wrf_oper.wrf_met(ruta_raiz=ruta_main,
                           ruta_metx=ruta_wpsx,
                           fecha_ini=fecha_ini,
                           fecha_fin=fecha_fin)
metgrid.copiar_archivos()
metgrid.editar_namelist()
metgrid.ejecutar_metgrid()

real = wrf_oper.wrf_real(ruta_raiz=ruta_main,
                         ruta_realx=ruta_arwx,
                         fecha_ini=fecha_ini,
                         fecha_fin=fecha_fin)
real.copiar_archivos()
real.editar_namelist()
real.ejecutar_real()

wrf = wrf_oper.wrf_wrf(ruta_raiz=ruta_main,
                       ruta_wrfx=ruta_arwx,
                       ruta_mpix=ruta_mpir,
                       nproc='64',
                       name_nodo='NODO-B-003',
                       fecha_ini=fecha_ini)
wrf.ejecutar_wrf()

copiar = wrf_oper.oper_copiar(ruta_raiz=ruta_main,
                              fecha_ini=fecha_ini)

copiar.copia_wrfout()
copiar.copia_wrfhr()
copiar.copia_wrfrst()
copiar.copia_config()
copiar.copia_gfs()
