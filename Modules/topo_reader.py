

def topo_file_id(topo_file):
    '''
    topo_type is a dictionary that standardises different variable
    naming of lon, lat, topo.
    Feel free to add your topography
    '''
    topo_type = {}
    topo_type['srtm']=None
    # etopo5.nc from Jeroen's easygrid
    if 'etopo5' in topo_file.split('/')[-1].lower():
        topo_type['lon'] = 'topo_lon'
        topo_type['lat'] = 'topo_lat'
        topo_type['topo'] = 'topo'
    # etopo2.nc from Romstools
    elif 'etopo2' in topo_file.split('/')[-1].lower():
        topo_type['lon'] = 'longitude'
        topo_type['lat'] = 'latitude'
        topo_type['topo'] = 'z'
    elif 'etopo1' in topo_file.split('/')[-1].lower():
        topo_type['lon'] = 'x'
        topo_type['lat'] = 'y'
        topo_type['topo'] = 'z'
    # srtm file
    elif 'srtm30' in topo_file.split('/')[-2].lower():
        topo_type['lon'] = 'longitude'
        topo_type['lat'] = 'latitude'
        topo_type['topo'] = 'elevation'
        topo_type['srtm']= True
    # croco file
    elif 'croco' in topo_file.split('/')[-1].lower():
        topo_type['lon'] = 'lon_rho'
        topo_type['lat'] = 'lat_rho'
        topo_type['topo'] = ''
    # gebco file
    elif 'gebco' in topo_file.split('/')[-1].lower():
        try: # unprocessed gebco file
            #nc = netcdf.Dataset(topo_file)
            #GEBCO in nc.title
            #nc.close()
            topo_type['lon']  = 'x_range'
            topo_type['lat']  = 'y_range'
            topo_type['topo'] = 'z'
        except:
            # processed gebco file
            topo_type['lon'] = 'lon'
            topo_type['lat'] = 'lat'
            topo_type['topo'] = 'topo'

    return topo_type

