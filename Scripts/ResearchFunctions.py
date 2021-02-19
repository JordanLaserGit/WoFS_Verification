import pandas
import os
from netCDF4 import Dataset
import numpy as np
import geopy.distance

#Inputs are the observation time and the path to the realstorm file.
def Stormloc(obstime,stormpath):

	print('Observation time: {}'.format(obstime))
	print('Path of real storm file: {}'.format(stormpath))

	storm = pandas.read_csv(stormpath)
	stormtime = np.asarray(storm.values[:,0])
	stormlat = np.asarray(storm.values[:,1], dtype = np.float32)
	stormlon = np.asarray(storm.values[:,2], dtype = np.float32)

	#if there is no time that matches the obstime time
    #initialize list for time distances
	timedistlist = []

    #time distance loop, finds how far obstime is from times in realstorm file
	for i in range(len(stormtime)):
	    timedist = abs(float(obstime) - stormtime[i])
	    timedistlist.append(timedist)

	#Locate closest time in storm file and use the lat/lon from this scan

	stormindex = np.where(timedistlist == np.nanmin(timedistlist))[0]

	diff = stormtime[stormindex] - float(obstime)
	print("Index is: " + str(stormindex))
	print("Difference is: " + str(diff))

	#Interpolating lat/lon if obstime is between scans, selecting scan time if obstime is a scan time
	if len(diff) == 1:
		if diff < 0:
		    Rlat = (stormlat[stormindex]*((5-abs(diff))/5) + (stormlat[stormindex+1]*((abs(diff))/5)))
		    Rlon = (stormlon[stormindex]*((5-abs(diff))/5) + (stormlon[stormindex+1]*((abs(diff))/5)))
		    print("obstime not found in storm time, interpolating")
		if diff > 0:
		    Rlat = (stormlat[stormindex-1]*((diff)/5) + (stormlat[stormindex]*(5-diff))/5)
		    Rlon = (stormlon[stormindex-1]*((diff)/5) + (stormlon[stormindex]*(5-diff))/5)
		    print("obstime not found in storm time, interpolating")
		if diff == 0:
			Rlat = stormlat[stormindex]
			Rlon = stormlon[stormindex]
			print("obstime found in storm time")

	# Averaging scan times lat/lon if obstime is exactly between two scan times
	if len(diff) == 2:
	    print(stormlat[stormindex])
	    Rlat = np.nanmean(stormlat[stormindex])
	    Rlon = np.nanmean(stormlon[stormindex])
	    print("obstime not found in storm time, interpolating")

	# Grab the first scan time if obstime is before all scan times, shouldn't happen
	if (float(obstime) <= stormtime[0]) & (float(obstime) >= 2000):
		Rlat = stormlat[0]
		Rlon = stormlon[0]
		print('Check scan times. obstime is before earliest scan time.')


	print("Storm lat: " + str(Rlat))
	print("Storm lon: " + str(Rlon))
	return Rlat, Rlon

# Inputs are storm lat/lon and obs lat/lon, calculates zonal and meridional displacement
# Requires geopy module
def RelativeCoordinates(stormlat,stormlon,obslat,obslon):

	coords_1 = (stormlat, stormlon)
	coords_2 = (obslat, obslon)

	xlat = float(obslat) - float(stormlat)
	xlon = float(obslon) - float(stormlon) 

	Re = 6.4*10**6
	latref = 40
	lonref = 100
	latconvert = 2*np.pi*Re/360
	lonconvert = 2*np.pi*Re*np.cos((np.pi/180)*latref)/360

	deltaY = (latconvert*xlat)/1000
	deltaX = (lonconvert*xlon)/1000

	distance = geopy.distance.distance(coords_1, coords_2).km

	return deltaX, deltaY, distance, xlat, xlon

# Inputs are the high resolution data, the altitude at which data exists, 
# and the altitudes that the high resolution data is decimated onto
# Returns the decimated data set
	
# Uses height as the vertical coordinate (used mostly for lidar, radiosonde kinematics)
def DecimateH(highres, highresalt, lowresalt):

	window = 50
	X = len(lowresalt)
	highresdecimated = np.ones(X)*np.nan
	for i in range(X):
		z = np.where((highresalt >= (lowresalt[i]-window)) & (highresalt <= (lowresalt[i]+window)))[0]
		highresdecimated[i] = np.nanmean(highres[z])

	return highresdecimated	
	
# Uses pressure as a vertical coordinate (used for radiosondes thermodynamics)
# This function can likely be simplified to resemble DecimateH, should work fine the way it is.
# I coded it this way to mitigate a pressure difference problem that was later resolved.
def DecimateP(highres, highresalt, lowresalt, highrespressure, lowrespressure):

	window = 50
	X = len(lowrespressure)
	highresdecimated = np.ones(X)*np.nan
	for i in range(X):
		z = []
		k = np.argmin(abs(lowrespressure[i] - highrespressure))
		diff = lowrespressure[i] - highrespressure[k]

# This code requires the pressure levels are within 5 hPa, otherwise highresdecimated remains nan for vertical level i
		if abs(diff) < 5:
			ind = np.argmin(abs(lowrespressure[i] - highrespressure))
			z = np.where((highresalt >= (highresalt[ind]-window)) & (highresalt <= (highresalt[ind]+window)))[0]
			highresdecimated[i] = np.nanmean(highres[z])
			
	return highresdecimated
