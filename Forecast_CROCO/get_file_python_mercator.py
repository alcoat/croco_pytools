######################################################################
#
#  Extract a subgrid from ECCO to get a CROCO forcing
#   Store that into monthly files.
#   Take care of the Greenwitch Meridian.
#
######################################################################
# Usage: 
#=======
# get_file_python_mercator(pathMotu,mercator_type,['varname1' 'varname2'],
#                          [lonmin lonmax latmin latmax, depthmin depthmax],
#                          ['startdate' 'enddate'],[user password],
#                          'outname.nc')
# 
# Check http://marine.copernicus.eu/web/34-products-and-services-faq.php 
# for changes in the command line: addresses, file names, variable names ...
#
# Currently needs motu-client.py v.1.0.8 and Python 2.7.x
import subprocess
from croco_tools_params import *

def get_mercator(pathmotu,mercator_type,url,sid,pid,vars,geom,date,info,outname):
	if mercator_type==1: # Mercator data 1/12 deg
		command = ['python -m motuclient --motu '+url+' --service-id '+sid+' --product-id '+pid+
		' --longitude-min '+geom[0]+' --longitude-max '+geom[1]+' --latitude-min '+geom[2]+' --latitude-max '+geom[3]+
		' --date-min '+date[0]+' --date-max '+date[1]+' --depth-min '+geom[4]+' --depth-max '+geom[5]+
		vars+' --out-dir ./ --out-name '+outname+' --user '+info[0]+' --pwd '+info[1]]
		print(command)
		p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

	else:              # UK Met Office data 1/4 deg
		command = 'python -m motuclient --motu '+\
			url+' --service-id '+sid+' --product-id '+pid+\
				' --longitude-min '+geom[0]+' --longitude-max '+geom[1]+' --latitude-min '+geom[2]+' --latitude-max '+geom[3]+\
					' --date-min '+date[0]+' --date-max '+date[1]+\
						' --depth-min '+geom[4]+' --depth-max '+geom[5]+\
							vars+'--out-dir ./ '+'--out-name '+outname+\
								' --user '+info[0]+' --pwd '+info[1]
		print(command)
		p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
	
	print('Waiting Download')
	(output, err) = p.communicate()
	print(output)
	p_status = p.wait()
	print('Download Finished')
	return


  
#===
# If you use proxy server, uncomment the following line and replace by your 
#   proxy url and port server. Beware that a SPACE is needed between 
#   "--proxy-server=" and "the proxy-server-name" !
#command[end+2]=sprintf('--proxy-server= http://your_proxy_server:your_proxy_port')
#===
	