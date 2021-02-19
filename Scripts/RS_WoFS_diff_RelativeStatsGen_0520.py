# Author: Jordan Laser
# This script calculates the differences between radiosonde observed variables and WoFS modeled variables
# Analysis is storm relative (requires storm locations, CSVs)
# WoFS and observed u,v,h,p,t,td,q are written to a netCDF
# Note: I use "decimating" throughout this script, perhaps "averaging" would have been a better word to use...

# Pat! Put the path to 'Pat_20th_rerun' here!
pat_path = 

# Import modules
import numpy as np
import netCDF4
from netCDF4 import Dataset
import pylab
import pandas as pd
import math as math
import glob
import os
import sys
from optparse import OptionParser

# Import scripts
from wofs_sounding_variable_extractor_RS_Relative_0517 import *
from ResearchFunctions import *

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

# Set the date
date = '0520'


############## read in observation files ################
files = os.listdir(indir)
files.sort()
print('Here are all the files in the folder you designated: ' + str(files))

for f, infile in enumerate(files):

	infile = os.path.join(indir, infile)
	df = pd.read_csv(infile, skiprows = 3, sep = ',', header = None)
	
	time = infile[-10:-6]

	
	print('Obs date: ' + str(date))
	print('Obs time: ' + str(time) + ' UTC')
	
	#Need example grid in order to determine closest surface grid point to observation
	Nfile = infile
	WoFSfile = sounding_extractor(time,date,Nfile)
	print('WoFS Sounding file used in comparison: ' + WoFSfile)
	

	#files date/times
	WoFSrun = WoFSfile[21:25]
	WoFSvalid = WoFSfile[26:30]
	print( WoFSrun + ' UTC' + ' forecast, valid at ' + WoFSvalid +  ' UTC')
	WoFS = Dataset('{}/Pat_20th_rerun/Radiosonde WoFS Sounding Files/{}'.format(pat_path, WoFSfile))

	#Pull wind components
	#WoFS variables have dimensions [member, z-level, x, y]
	Sxdex, Sydex = [1,1] #pulling middle grid point from the 9 provided in the sounding file
	Wlat = WoFS.variables['xlat']
	Wlon = WoFS.variables['xlon']
	uWoFS = WoFS.variables['u'][:,:,Sxdex,Sydex]
	vWoFS = WoFS.variables['v'][:,:,Sxdex,Sydex]
	tWoFS = WoFS.variables['t'][:,:,Sxdex,Sydex] +273.15
	pWoFS = WoFS.variables['p'][:,:,Sxdex,Sydex] 
	hWoFS = WoFS.variables['z_agl'][:,:,Sxdex,Sydex]
	qWoFS = WoFS.variables['q'][:,:,Sxdex,Sydex]
	
	thetaWoFS = (tWoFS+273.15)*(1000/pWoFS)**(287/1004)
	thetaWoFS = thetaWoFS - 273.15
	
	
	member_list = [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 2, 3, 4, 5, 6, 7, 8, 9]
	#ysu = [1,2,7,8,13,14]
	#myj = [3,4,9,10,15,16]            #pbl schemes based off of member
	#mynn = [5,6,11,12,17,18]
	ysu = [0,4,5,10,15,16]
	myj = [1,6,7,11,12,17]             #pbl indicies
	mynn = [2,3,8,9,13,14]

	#radiosonde variables
	#set columns of df to variable names
	temp = df[2]
	rh = df[3]
	dew = df[4]
	pressure = df[5]
	winddir = df[6]
	windspd = df[7]
	elv = df[8]
	alt = elv-elv[0]
	Slat = df[9]
	Slon = df[10]
	windU = -windspd*np.cos((winddir*3.1415/180)-(3.1415/2))
	windV = windspd*np.sin((winddir*3.1415/180)-(3.1415/2))
	w = df[13]
	
	#Calculate Potential Temperature (temp = kelvin)
	thetaObs = temp *(1000/pressure)**(287/1004)
	thetaObs = thetaObs - 273.15

	#Calculate water vapor mixing ratio from obs
	e_0 = 6.1173
	t_0 = 273.16
	Rv = 461.5
	Lv = 2501000
	es = e_0*np.exp((Lv/Rv)*((1/t_0)-(1/temp)))
	e = (rh/100)*es
	q = 0.622*(e/(pressure-e))	

	###################################################################
	#             DECIMATING OBSERVATION ONTO WOFS GRID
	###################################################################

	deciU = DecimateH(windU,alt,hWoFS[0,:])
	deciV = DecimateH(windV,alt,hWoFS[0,:])
	decit = DecimateP(temp,alt,hWoFS[0,:],pressure,pWoFS[0,:])
	deciq = DecimateP(q,alt,hWoFS[0,:],pressure,pWoFS[0,:])
	decidew = DecimateP(dew,alt,hWoFS[0,:],pressure,pWoFS[0,:])
	

	###################   STATISTICS   ########################
	#calculating difference between observation and each member, so diff should have dimensions of vertical level vs member
	#Kiematic Differences
	#Zonal Component
	diffU = np.ones((18,len(deciU)))*np.nan
	for i in range(18):
		for j in range(len(deciU)):
			diffU[i,j] = uWoFS[i,j]-deciU[j]

	#Meridional Component
	diffV = np.ones((18,len(deciV)))*np.nan
	for i in range(18):
		for j in range(len(deciV)):
			diffV[i,j] = vWoFS[i,j]-deciV[j]


	#Thermodynamic differences

	diffT = np.ones((18,len(decit)))*np.nan
	for i in range(18):
		for j in range(len(decit)):
			diffT[i,j] = (tWoFS[i,j]-decit[j])

	diffq = np.ones((18,len(deciq)))*np.nan
	for i in range(18):
		for j in range(len(deciq)):
			diffq[i,j] = qWoFS[i,j]-deciq[j]


	###########################################################################
	###########################################################################
	############################## Writing output file ########################
	ne = 18
	nz = len(deciU)

	outdir = "{}/Pat_20th_rerun/Radiosonde Stat Files/".format(pat_path)
	outname = "WoFSMemberStats_{}_{}_{}.nc".format(date,WoFSrun,time)         #output file
	output_path = os.path.join(outdir,outname)
	try:
		fout = netCDF4.Dataset(output_path, "w")
	except:
	   print("Could not create %s!\n" % output_path)


	### Create file and dimensions: ###
	fout.createDimension('NE', ne)
	fout.createDimension('NZ', nz)

	### Create variables ###
	diffu_var = fout.createVariable('DU', 'f4', ('NE','NZ',))
	diffu_var.long_name = "difference in U-component of Wind"
	diffu_var.units = "m s**-1"

	diffv_var = fout.createVariable('DV', 'f4', ('NE','NZ',))
	diffv_var.long_name = "Difference in V-component of Wind"
	diffv_var.units = "m s**-1"

	difft_var = fout.createVariable('DT', 'f4', ('NE','NZ',))
	difft_var.long_name = "Difference in Temperature"
	difft_var.units = "\u00B0C"

	diffq_var = fout.createVariable('Dq', 'f4', ('NE','NZ',))
	diffq_var.long_name = "Difference in Mixing Ratio"
	diffq_var.units = "Kg/Kg"
	
	
	
	uwofs_var = fout.createVariable('uwofs', 'f4', ('NE','NZ',))
	uwofs_var.long_name = "WoFS U wind"
	uwofs_var.units = "m/s"
	
	vwofs_var = fout.createVariable('vwofs', 'f4', ('NE','NZ',))
	vwofs_var.long_name = "WoFS V wind"
	vwofs_var.units = "m/s"

	thetaWoFS_var = fout.createVariable('thetaWoFS', 'f4', ('NE','NZ',))
	thetaWoFS_var.long_name = "WoFS Temperature"
	thetaWoFS_var.units = "K"
	
	tWoFS_var = fout.createVariable('tWoFS', 'f4', ('NE','NZ',))
	tWoFS_var.long_name = "WoFS Temperature"
	tWoFS_var.units = "K"
	
	qWoFS_var = fout.createVariable('qWoFS', 'f4', ('NE','NZ',))
	qWoFS_var.long_name = "WoFS Water Vapor Mixing Ratio"
	qWoFS_var.units = "kg/kg"
	
	pressureWoFS_var = fout.createVariable('pWoFS', 'f4', ('NE','NZ',))
	pressureWoFS_var.long_name = "WoFS Pressure"
	pressureWoFS_var.units = "hPa"
	
	diffq_var = fout.createVariable('h', 'f4', ('NE','NZ',))
	diffq_var.long_name = "Height of Veritcal Levels in WoFS"
	diffq_var.units = "m"
	
	
	
	deciU_var = fout.createVariable('deciU', 'f4', ('NZ',))
	deciU_var.long_name = "Decimated Observed U-component of Wind"
	deciU_var.units = "m/s"

	deciV_var = fout.createVariable('deciV', 'f4', ('NZ',))
	deciV_var.long_name = "Decimated Observed V-component of Wind"
	deciV_var.units = "m/s"

	decit_var = fout.createVariable('decit', 'f4', ('NZ',))
	decit_var.long_name = "Decimated Observed Temperature"
	decit_var.units = "K"

	deciq_var = fout.createVariable('deciq', 'f4', ('NZ',))
	deciq_var.long_name = "Decimated Observed Mixing Ratio"
	deciq_var.units = "Kg/Kg"
	
	dew_var = fout.createVariable('Td', 'f4', ('NZ',))
	dew_var.long_name = "Observed Dewpoint Temperature"
	dew_var.units = "K"
	
	
	###Write Variables###
	fout.variables['DU'][:] = diffU
	fout.variables['DV'][:] = diffV
	fout.variables['DT'][:] = diffT
	fout.variables['Dq'][:] = diffq
	
	fout.variables['uwofs'][:] = uWoFS
	fout.variables['vwofs'][:] = vWoFS
	fout.variables['thetaWoFS'][:] = thetaWoFS
	fout.variables['tWoFS'][:] = tWoFS
	fout.variables['qWoFS'][:] = qWoFS
	fout.variables['pWoFS'][:] = pWoFS
	fout.variables['h'][:] = hWoFS
	
	fout.variables['deciU'][:] = deciU
	fout.variables['deciV'][:] = deciV
	fout.variables['decit'][:] = decit
	fout.variables['deciq'][:] = deciq
	fout.variables['Td'][:] = decidew
	

	fout.close()
	del fout