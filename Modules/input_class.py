import numpy as np
import tools
import netCDF4 as netcdf
import inputs_readers as dico
'''
This class need to:
    - open netcdf
    - read depth,lon,lat(T,U,V if needed),ssh, temp, salt, u,v
'''

class getdata():   
    def __init__(self,inputdata,inputfile,crocogrd,bdy=None): # bdy=[obs,tstart,tend,cycle]
        
        self.var=dico.lookvar(inputdata) # Dictionary to find the names of the input variables
        if bdy is None: # Ini case
            if inputdata != 'eccov4':
                self.depth=tools.read_nc(inputfile,self.var['depth'])
  
                self.ncglo   = { 'ssh'  : netcdf.Dataset(inputfile,'r').variables,\
                                 'temp' : netcdf.Dataset(inputfile,'r').variables,\
                                 'salt' : netcdf.Dataset(inputfile,'r').variables,\
                                 'u'    : netcdf.Dataset(inputfile,'r').variables,\
                                 'v'    : netcdf.Dataset(inputfile,'r').variables\
                               }
        
            else:
                self.depth=tools.read_nc(inputfile['temp'],self.var['depth']) # read depth in temp file (they are the same in all files)
            
                self.ncglo   = { 'ssh'  : netcdf.Dataset(inputfile['ssh'],'r').variables,\
                                 'temp' : netcdf.Dataset(inputfile['temp'],'r').variables,\
                                 'salt' : netcdf.Dataset(inputfile['salt'],'r').variables,\
                                 'u'    : netcdf.Dataset(inputfile['u'],'r').variables,\
                                 'v'    : netcdf.Dataset(inputfile['v'],'r').variables\
                               }

            [self.lonT ,self.latT ,self.idmin  ,self.idmax  ,self.jdmin  ,self.jdmax  ,self.period  ]  = self.handle_periodicity(crocogrd,'r')
            [self.lonU ,self.latU ,self.idminU ,self.idmaxU ,self.jdminU ,self.jdmaxU ,self.periodU ]  = self.handle_periodicity(crocogrd,'u')
            [self.lonV ,self.latV ,self.idminV ,self.idmaxV ,self.jdminV ,self.jdmaxV ,self.periodV ]  = self.handle_periodicity(crocogrd,'v')        
        elif bdy is not None and bdy[-1]==0: # bdy case
            self.depth=tools.read_nc(inputfile[0],self.var['depth'])
            
            if inputdata != 'eccov4':
                self.ncglo   = { 'ssh'  : netcdf.MFDataset(inputfile,'r',aggdim=self.var['time_dim']).variables,\
                                 'temp' : netcdf.MFDataset(inputfile,'r',aggdim=self.var['time_dim']).variables,\
                                 'salt' : netcdf.MFDataset(inputfile,'r',aggdim=self.var['time_dim']).variables,\
                                 'u'    : netcdf.MFDataset(inputfile,'r',aggdim=self.var['time_dim']).variables,\
                                 'v'    : netcdf.MFDataset(inputfile,'r',aggdim=self.var['time_dim']).variables,\
                                 'time' : netcdf.MFDataset(inputfile,'r',aggdim=self.var['time_dim']).variables\
                               }
            else:
                self.depth=netcdf.MFDataset(inputfile['temp'],'r',aggdim=self.var['time_dim']).variables[self.var['depth']][:] # read depth in temp file (they are the same in all files)

                self.ncglo   = { 'ssh'  : netcdf.MFDataset(inputfile['ssh'],'r',aggdim=self.var['time_dim']).variables,\
                                 'temp' : netcdf.MFDataset(inputfile['temp'],'r',aggdim=self.var['time_dim']).variables,\
                                 'salt' : netcdf.MFDataset(inputfile['salt'],'r',aggdim=self.var['time_dim']).variables,\
                                 'u'    : netcdf.MFDataset(inputfile['u'],'r',aggdim=self.var['time_dim']).variables,\
                                 'v'    : netcdf.MFDataset(inputfile['v'],'r',aggdim=self.var['time_dim']).variables,\
                                 'time' : netcdf.MFDataset(inputfile['ssh'],'r',aggdim=self.var['time_dim']).variables\
                               }
   

            for boundary, is_open in zip(bdy[0].keys(), bdy[0].values()):

                if 'west' in boundary and is_open:
                    print('\nHandling western grid')
                    print('---------------------')
                    [self.lonTW ,self.latTW ,self.idminW ,self.idmaxW ,self.jdminW ,self.jdmaxW ,self.periodW  ]  = self.handle_periodicity(crocogrd,'r',bdy='west')
                    [self.lonUW ,self.latUW ,self.idminUW ,self.idmaxUW ,self.jdminUW ,self.jdmaxUW ,self.periodUW ]  = self.handle_periodicity(crocogrd,'u',bdy='west')
                    [self.lonVW ,self.latVW ,self.idminVW ,self.idmaxVW ,self.jdminVW ,self.jdmaxVW ,self.periodVW ]  = self.handle_periodicity(crocogrd,'v',bdy='west')
                elif 'east' in boundary and is_open:
                    print('\nHandling eastern grid')
                    print('---------------------')
                    [self.lonTE ,self.latTE ,self.idminE ,self.idmaxE ,self.jdminE ,self.jdmaxE ,self.periodE  ]  = self.handle_periodicity(crocogrd,'r',bdy='east')
                    [self.lonUE ,self.latUE ,self.idminUE ,self.idmaxUE ,self.jdminUE ,self.jdmaxUE ,self.periodUE ]  = self.handle_periodicity(crocogrd,'u',bdy='east')
                    [self.lonVE ,self.latVE ,self.idminVE ,self.idmaxVE ,self.jdminVE ,self.jdmaxVE ,self.periodVE ]  = self.handle_periodicity(crocogrd,'v',bdy='east')
                elif 'south' in boundary and is_open:
                    print('\nHandling southern grid')
                    print('----------------------')
                    [self.lonTS ,self.latTS ,self.idminS ,self.idmaxS ,self.jdminS ,self.jdmaxS ,self.periodS  ]  = self.handle_periodicity(crocogrd,'r',bdy='south')
                    [self.lonUS ,self.latUS ,self.idminUS ,self.idmaxUS ,self.jdminUS ,self.jdmaxUS ,self.periodUS ]  = self.handle_periodicity(crocogrd,'u',bdy='south')
                    [self.lonVS ,self.latVS ,self.idminVS ,self.idmaxVS ,self.jdminVS ,self.jdmaxVS ,self.periodVS ]  = self.handle_periodicity(crocogrd,'v',bdy='south')
                elif 'north' in boundary and is_open:
                    print('\nHandling northern grid')
                    print('----------------------')
                    [self.lonTN ,self.latTN ,self.idminN ,self.idmaxN ,self.jdminN ,self.jdmaxN ,self.periodN  ]  = self.handle_periodicity(crocogrd,'r',bdy='north')
                    [self.lonUN ,self.latUN ,self.idminUN ,self.idmaxUN ,self.jdminUN ,self.jdmaxUN ,self.periodUN ]  = self.handle_periodicity(crocogrd,'u',bdy='north')
                    [self.lonVN ,self.latVN ,self.idminVN ,self.idmaxVN ,self.jdminVN ,self.jdmaxVN ,self.periodVN ]  = self.handle_periodicity(crocogrd,'v',bdy='north')

              
    #####################################
    def handle_periodicity(self,crocogrd,grid,bdy=None):
        '''
        handle_periodicity checks whether domain is inside the croco file.
        If so, check if there is a need for create a periodicity between
        the last and first longitude points ( for global data).
        It is returning lon/lat/topo adapted to the desired domain
        geolim = [lonmin,lonmax,latmin,latmax]
 
        input: inputfile : data file
               crocogrd  : croco grid lon/lat
               grid      : which grid (rho: 'r', u: 'u', v: 'v')
 
        output: lon/lat of the grid
                imin/imax index min/max xaxis
                jmin/jmax index min/max yaxis
        '''

        print('Reading coordinate file for grid: '+grid )
        print('-----------------------------------')
        if grid == 'r':
            inpgrd = 'ssh'
        else:
            inpgrd = grid
        lon=self.ncglo[inpgrd][self.var['lon'+grid]]
        lat=self.ncglo[inpgrd][self.var['lat'+grid]]

        if len(lon.shape)==2: # Some datasets are a bit different
            lon = self.ncglo[inpgrd][self.var['lon'+grid]][0,:]
            lat = self.ncglo[inpgrd][self.var['lat'+grid]][:,0]

        for i in range(1,lon.shape[0]): # Fix discontinuity
            if lon[i]<lon[i-1]:        # between 180/-180 in the
                lon[i]=lon[i]+360      # middle
        ####
        if bdy is not None:
            geolim=[np.min(eval(''.join(('crocogrd.lon_',bdy))))-1,np.max(eval(''.join(('crocogrd.lon_',bdy))))+1,\
                    np.min(eval(''.join(('crocogrd.lat_',bdy))))-1,np.max(eval(''.join(('crocogrd.lat_',bdy))))+1]
        else:
            geolim=[crocogrd.lonmin()-1,crocogrd.lonmax()+1,crocogrd.latmin()-1,crocogrd.latmax()+1]

        jmin=tools.indx_bound(lat, geolim[2])
        jmax=tools.indx_bound(lat, geolim[-1])

        if 0 < jmin and jmin < lat.shape[0] and 0 < jmax and jmax < lat.shape[0] :
            if jmin > 1 :
                jmin=jmin-1
            jmax=jmax+2
        else:
            print('North-south extents of the dataset ',lat[0],lat[-1],' are not sufficient to cover the entire model grid.')
            exit()
        ####
        imin=tools.indx_bound(lon, geolim[0])
        imax=tools.indx_bound(lon, geolim[1])

        if 0 < imin and imin < lon.shape[0] and 0 < imax and imax < lon.shape[0] :
            if imax > 1:
                imin=imin-1
            imax=imax+2
            shft_west=0 ; shft_east=0 ; period=0
            print('Single region dataset imin/imax=',imin,imax, )
        else:
        ######
            ptest=lon[-1]-lon[0]-360
            dx=(lon[-1]-lon[0])/(lon.shape[0]-1)
            epsil=0.01*abs(dx)
            if abs(ptest) < epsil :
                period=lon.shape[0]-1
            elif abs(ptest+dx) < epsil :
                period=lon.shape[0]
            else:
                period=0

            if period>0:
                print('Identified periodicity domain in data of ', period,' points out of', lon.shape[0])
            else :
                print('ERROR: The data does not cover the entire grid. Change your grid definition')
                exit()
        ##
            shft_west=0
            if imin==0 :
                shft_west=-1
                imin=tools.indx_bound(lon, geolim[0]+360)
            elif imin==lon.shape[0] :
                shft_west=+1
                imin=tools.indx_bound(lon, geolim[0]-360)
        ##
            shft_east=0
            if imax == 0:
                shft_east=-1
                imax=tools.indx_bound(lon, geolim[1]+360)
            elif imax == lon.shape[0]:
                shft_east=+1
                imax=tools.indx_bound(lon, geolim[1]-360)
    
            if 0<imin and imin <lon.shape[0] and 0<imax and imax<lon.shape[0] :
                if imin>1:
                    imin=imin-1
                imax=imax+1
            else:
                print('ERROR: Data longitude covers 360 degrees, but still cannot find  starting and ending indices.')
                exit()

        print('Bounding indices of the relevant part to be extracted from the entire dataset:\n', 
              'imin,imax =', imin,imax,'out of', lon.shape[0],'jmin,jmax =',jmin,jmax, 'out of',lat.shape[0])
        ny_lat=jmax-jmin+1
        start2=jmin ; end2=start2+ny_lat; count2=ny_lat
        lat_tmp=np.zeros([ny_lat])
        for j in range(0,ny_lat):
            lat_tmp[j]=lat[j+jmin-1]
        #####
        if imin < imax :
            nx_lon=imax-imin+1
            start1=imin ; end1=start1+nx_lon ; count1=nx_lon

            ishft=imin-1
            lon_tmp=np.zeros([nx_lon])
            if shft_west>0 and shft_east>0:
                for i in range(0,nx_lon):
                    lon_tmp[i]=lon[i+ishft] +360
            elif shft_west<0 and shft_east<0:
                for i in range(0,nx_lon):
                     lon_tmp[i]=lon[i+ishft]-360
            elif shft_west== 0 and shft_east==0:
                for i in range(0,nx_lon) :
                    lon_tmp[i]=lon[i+ishft]
            else:
                print('Error in shifting algoritm')
                exit()
            (lon,lat)=np.meshgrid(lon_tmp,lat_tmp)
        ###
        elif imin>imax:
            print('Reading topography in two separate parts adjacent through 360-degree periodicity\n First...' )
            nx_lon=imax+period-imin+1
            xtmp  = np.zeros([nx_lon])
            start1=0 ; end1=start1+nx_lon; count1=imax

            ishft=nx_lon-count1
            if shft_east>0:
                for i in range(0,count1):
                    xtmp[i+ishft]=lon[i] +360
            elif shft_east<0:
                for i in range(0,count1):
                    xtmp[i+ishft]=lon[i] -360
            else:
                for i in range(0,count1):
                    xtmp[i+ishft]=lon[i]

            print('Second...')
            start1=imin ; count1=period-imin; end1=start1+count1
            ishft=imin-1
            if shft_west>0:
                for i in range(0,count1):
                    xtmp[i]=lon[i+ishft] +360
            elif shft_west<0 :
                for i in range(0,count1):
                    xtmp[i]=lon[i+ishft] -360
            else:
                for i in range(0,count1):
                    xtmp[i]=lon[i+ishft]
            lon_tmp=np.zeros([xtmp.shape[0]])
            for i in range(0,nx_lon):
                lon_tmp[i]=xtmp[i]

            del lon,lat
            (lon,lat)=np.meshgrid(lon_tmp,lat_tmp)
    
        return lon,lat,imin,imax,jmin,jmax,period

    #############################   
    def var_periodicity(self,vname,l,k,bdy=""):
        '''
        handle periodicity for tracers. Limits (imin,imax,jmin,jmax) are fixed before
        by runing 
        '''
        if vname == 'u':
            imin=eval(''.join(("self.idminU"+bdy))) ; imax=eval(''.join(("self.idmaxU"+bdy)))
            jmin=eval(''.join(("self.jdminU"+bdy))) ; jmax=eval(''.join(("self.jdmaxU"+bdy)))
            period=eval(''.join(("self.periodU"+bdy)))
        elif vname == 'v':
            imin=eval(''.join(("self.idminV"+bdy))) ; imax=eval(''.join(("self.idmaxV"+bdy)))
            jmin=eval(''.join(("self.jdminV"+bdy))) ; jmax=eval(''.join(("self.jdmaxV"+bdy)))
            period=eval(''.join(("self.periodV"+bdy)))
        else:
            imin=eval(''.join(("self.idmin"+bdy))) ; imax=eval(''.join(("self.idmax"+bdy)))
            jmin=eval(''.join(("self.jdmin"+bdy))) ; jmax=eval(''.join(("self.jdmax"+bdy)))
            period=eval(''.join(("self.period"+bdy)))

        ny_lat=jmax-jmin+1
        start2=jmin ; end2=start2+ny_lat; count2=ny_lat

        if imin < imax :
            nx_lon=imax-imin+1
            start1=imin ; end1=start1+nx_lon ; count1=nx_lon

            if k==-1:
                field=np.array(self.ncglo[vname][self.var[vname]][l,start2:end2,start1:end1])
            else:
                field=np.array(self.ncglo[vname][self.var[vname]][l,k,start2:end2,start1:end1])

        elif imin>imax:    
            nx_lon=imax+period-imin+1
            try:
                lent=l.shape[0]
                ftmp = np.zeros([lent,ny_lat,nx_lon])
            except:
                ftmp = np.zeros([ny_lat,nx_lon])
            # First
            start1=0 ; end1=start1+nx_lon; count1=imax
        
            if k==-1:
                field=np.array(self.ncglo[vname][self.var[vname]][l,start2:end2,start1:end1])
            else:
                field=np.array(self.ncglo[vname][self.var[vname]][l,k,start2:end2,start1:end1])

            if len(ftmp.shape)<3:
                for j in range(0,count2):
                    for i in range(0,count1):
                        ftmp[j,nx_lon-imax+i-1]=field[j,i]
            else:
                for j in range(0,count2):
                    for i in range(0,count1):
                        ftmp[:,j,nx_lon-imax+i-1]=field[:,j,i]

            del field

            # Second
            start1=imin ; count1=period-imin; end1=start1+count1
            if k==-1:
                field=np.array(self.ncglo[vname][self.var[vname]][l,start2:end2,start1:end1])
            else:
                field=np.array(self.ncglo[vname][self.var[vname]][l,k,start2:end2,start1:end1])

            if len(ftmp.shape)<3:
                for j in range(0,count2):
                    for i in range(0,count1):
                        ftmp[j,i]=field[j,i]
            else:
                for j in range(0,count2):
                    for i in range(0,count1):
                        ftmp[:,j,nx_lon-imax+i-1]=field[:,j,i]

            del field

            field=np.copy(ftmp)

        return field



