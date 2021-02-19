import numpy as np
import netCDF4
from netCDF4 import Dataset
import numpy.ma as ma
import pandas as pd
import math as math
import glob
import xarray
import os
import sys
from optparse import OptionParser
from wofs_sounding_variable_extractor_Lidar_0520 import *
from ResearchFunctions import *

# Pat! Put the path to 'Pat_20th_rerun' here!
pat_path = 

#################### Command Options ####################
parser = OptionParser()
parser.add_option("-o", dest="indir", type="string", default= None, help="Input file of observation to be compared with WoFS")

(options, args) = parser.parse_args()

if (options.indir == None):
   parser.print_help()
   print()
   sys.exit(1)
else:
   indir = options.indir


############## read in observation files ################
files = os.listdir(indir)
files.sort()
print('Here are all the files in the folder you designated: ' + str(files))

for f, infile in enumerate(files):
	infile = os.path.join(indir, infile)
	Nfile = infile
	Lidar = Dataset(infile)
	date = '0520'
	print('Obs date: ' + date)
	
	#Converts Lidar time output to HHMM format
	tLid = Lidar.variables['hour']
	tUTC = np.ones(len(tLid), dtype = np.int)
	minUTC = np.ones(len(tLid), dtype = np.int)
	for i in range(len(tLid)):
		minUTC[i] = int((tLid[i]-int(tLid[i]))*60)
		if minUTC[i] < 10:
			tUTC[i] = str(int(tLid[i])) + '0' + str(minUTC[i])
		else:
			tUTC[i] = str(int(tLid[i])) + str(minUTC[i])

	######## Deployment 1 and 2
	time = [2030,2035,2040,2045,2050,2055,2100,2105,2110,2115,2120,2125,2130,2135,2155,2200,2205,2210]

	#Need example grid in order to determine closest surface grid point to observation
	for x in range(len(time)):
		print(time[x])
		
		WoFSfile, hour = sounding_extractor(str(time[x]),date,Nfile)

		print('WoFS Sounding file used in comparison: ' + WoFSfile)
	
		#files date/time
		WoFSrun = WoFSfile[21:25]
		WoFSvalid = WoFSfile[26:30]
		print( WoFSrun + 'UTC' + ' forecast, valid at ' + WoFSvalid +  'UTC')
		WoFS = Dataset('{}/Pat_20th_rerun/Lidar WoFS Sounding Files/'.format(pat_path,WoFSfile))

		#WoFS variables have dimensions [member, z-level, x, y]
		xdex, ydex = [1,1] #pulling middle grid point from the 9 provided in the sounding file
		Wlat = WoFS.variables['xlat']
		Wlon = WoFS.variables['xlon']
		uWoFS = WoFS.variables['u'][:,:,xdex,ydex]
		vWoFS = WoFS.variables['v'][:,:,xdex,ydex] 
		hWoFS = WoFS.variables['z_agl'][:,:,xdex,ydex]
		member_list = [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 2, 3, 4, 5, 6, 7, 8, 9]
	
		#Lidar variables
		#Lidar wind variable have dimensions [t,z]
		uLid = Lidar.variables['u']
		vLid = Lidar.variables['v']
		hLid = Lidar.variables['height']
	
			
		# #Mean WoFS U,V,h
		uMeanWoFS = np.ones(len(uWoFS[1,:]))*np.nan
		for i in range(len(uWoFS[1,:])):
			uMeanWoFS[i] = np.average(uWoFS[:,i])
		vMeanWoFS = np.ones(len(vWoFS[1,:]))*np.nan
		for i in range(len(vWoFS[1,:])):
			vMeanWoFS[i] = np.average(vWoFS[:,i])
		hMeanWoFS = np.ones(len(hWoFS[1,:]))*np.nan
		for i in range(len(hWoFS[1,:])):
			hMeanWoFS[i] = np.average(hWoFS[:,i])
		
		#Averaging Lidar data based off of time
		Index = np.where((tUTC <= (int(time[x]) + 3)) & (tUTC >= (int(time[x]) - 3)))[0]
		print("Scan times used in averaging: " + str(tUTC[Index]))
		
		Uavg = np.nanmean(uLid[Index,:], axis = 0)
		Vavg = np.nanmean(vLid[Index,:], axis = 0)
		

	#########################################################
		# Decimating Lidar onto WoFS vertical grid
		
		deciU = DecimateH(Uavg,hLid,hWoFS[0,:])
		deciV = DecimateH(Vavg,hLid,hWoFS[0,:])


	###################   STATISTICS   #######################
		#calculating difference between observation and each member, 
		# so diff should have dimensions of vertical level vs member
		#Differences
		#Zonal Component
		diffU = np.ones((18,50))*np.nan
		for i in range(18):
			for j in range(50):
				diffU[i,j] = uWoFS[i,j]-deciU[j]

		#Meridional Component
		diffV = np.ones((18,50))*np.nan
		for i in range(18):
		
			for j in range(50):
				diffV[i,j] = vWoFS[i,j]-deciV[j]

	############# Writing .nc file #########################

		ne = 18
		nz = 50
		outdir = "{}/Pat_20th_rerun/Lidar Stat Files/".format(pat_path)
		outname = "Lidar_WoFSMemberStats_{}_{}_{}.nc" .format(date,hour,time[x])        #output file
		output_path = os.path.join(outdir,outname)
		try:
			fout = netCDF4.Dataset(output_path, "w")
		except:
		   print("Could not create %s!\n" % output_path)

		### Create file and dimensions: ###
		fout.createDimension('NE', ne)
		fout.createDimension('NZ', nz)

		### 3-D Variables ###
		### 3-D Variables ###
		deciu_var = fout.createVariable('DeciU', 'f4', ('NZ',))
		deciu_var.long_name = "decimated U-component of wind"
		deciu_var.units = "m s**-1"
		
		diffv_var = fout.createVariable('DeciV', 'f4', ('NZ',))
		diffv_var.long_name = "decimated V-component of wind"
		diffv_var.units = "m s**-1"
		
		uwofs_var = fout.createVariable('uWoFS', 'f4', ('NE','NZ',))
		uwofs_var.long_name = "WoFS U wind"
		uwofs_var.units = "m/s"
	
		vwofs_var = fout.createVariable('vWoFS', 'f4', ('NE','NZ',))
		vwofs_var.long_name = "WoFS V wind"
		vwofs_var.units = "m/s"
		
		diffu_var = fout.createVariable('DU', 'f4', ('NE','NZ',))
		diffu_var.long_name = "difference in U-component of wind"
		diffu_var.units = "m s**-1"

		diffv_var = fout.createVariable('DV', 'f4', ('NE','NZ',))
		diffv_var.long_name = "difference in V-component of wind"
		diffv_var.units = "m s**-1"
	
		
		h_var = fout.createVariable('h', 'f4', ('NE','NZ',))
		h_var.long_name = "Height of Veritcal Levels in WoFS"
		h_var.units = "m"
	
		###Write Variables###
		fout.variables['DeciU'][:] = deciU[:]
		fout.variables['DeciV'][:] = deciV[:]
		fout.variables['uWoFS'][:] = uWoFS
		fout.variables['vWoFS'][:] = vWoFS
		fout.variables['DU'][:] = diffU[:,:]
		fout.variables['DV'][:] = diffV[:,:]
		fout.variables['h'][:,:] = hWoFS[:,:]

		fout.close()
		del fout