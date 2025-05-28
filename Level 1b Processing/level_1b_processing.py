#!/usr/bin/env python
# coding: utf-8

# Project: Roman Space Telescope's Wide Field Instrument Integration and Test
# Author: Nargess Memarsadeghi, NASA GSFC, Code 586, Spring 2021.
# Description: This routine takes the complete path to a RST Level 1A .h5 input file and an output directory as command line inputs, and outputs level-1b data product as described in Roman Space Telescope's Wide Field Instrument Integration and Test Data Format Control Book.
#Example: $ python3 level_1b_processing.py L1A_input_filename output_directory mask_flag compression_flag
# sample file: /Users/nmemarsa/Documents/Projects/WFIRST/WFIRST_Data/L1A/20190916_95k_0p6m0p1_noise_20496_001.h5

import numpy as np
import sys, math, concurrent.futures, time #progress bar
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import h5py
import ntpath
import os

def readRSTFitsFile():
    global hdu, dataCube, dims, time_axis
    
    # get hdu information
    hdu = fits.open(inputFile)

    nreads = hdu[0].header.get('NREADS')
    print('nreads: ', nreads)

    # get the time_axis
    frtime = hdu[0].header.get('FRTIME')
    frskip = hdu[0].header.get('FRSKIP')
    cdstime = frtime * (frskip + 1)
    time_axis=(np.arange(nreads)+1)*cdstime
    print('time axis: ', time_axis)

    # Load image, subtract each frame from 2^16 -1 (except first one?), 
    print("Loading Images into memory...\n")
    dataCube = fits.getdata(inputFile)
    dims = dataCube.shape
    print('original data size using fits.getData: ', dims )
    print('\n');

    # get rid of extra singular dimension
    dataCube = np.squeeze(dataCube)

    # correct cube by substracting values from 2^16-1 or 65535
    dataCube= 65535.- dataCube.astype(np.float32)

    dataCube = np.array(dataCube)
    dims = dataCube.shape

    print('data size after squeeze operation: ', dims )
    print('dims[0]: ', dims[0])
    print('dims[1]: ', dims[1])
    print('dims[2]: ', dims[2])
    print('\n');

def readRST_L1A_File():
    global dataCube, dims, time_axis, root_attrs, frames_attrs, inFile
        
    #get root_attrs
    root_attrs = inFile.attrs
    #filedate = root_attrs['time_end_str'].decode()
    #print('*** filedate: ***', filedate)
    
    #root_attr_keys = list(inFile.attrs.keys())
    #root_attr_vals = list(inFile.attrs.values())
  
    
    dataCube = inFile.get('Frames')
    dims = dataCube.shape
    nreads = dims[0]
    print('nreads: ', nreads)
    
    #get frame_attrs
    frames_attrs = dataCube.attrs
    
    #frames_attr_keys = list(dataCube.attrs.keys())
    #frames_attr_vals = list(dataCube.attrs.values())
    
    # get the time_axis
    exptime = frames_attrs['exposure_length_actual']
    frtime = exptime /nreads
    frskip = 0;
    cdstime = frtime * (frskip + 1)
    time_axis=(np.arange(nreads)+1)*cdstime
    print('time axis: ', time_axis)

    # get rid of extra singular dimension
    dataCube = np.squeeze(dataCube)

    # No need for L1A data anymore to correct cube by substracting values from 2^16-1 or 65535
    # dataCube= 65535.- dataCube.astype(np.float32)

    dataCube = np.array(dataCube)
    dims = dataCube.shape

    print('data size after squeeze operation: ', dims )
    print('dims[0]: ', dims[0])
    print('dims[1]: ', dims[1])
    print('dims[2]: ', dims[2])
    print('\n');
    
  #  inFile.close()

# cds image calculation
def calc_cdsImage():
    global cdsImage
    cdsImage = np.squeeze(dataCube[dims[0]-1,:,:] - dataCube[0,:,:])

# Slope image calculation
# calculated from equation in https://en.wikipedia.org/wiki/Simple_linear_regression
def calc_raw_slopeImage():
    print('in calc_raw_slopeImage()\n')
    global slopeImage, time_axis
    
    time_axis = time_axis - np.mean(time_axis, axis=0)
    denominator = np.sum(time_axis**2)
    time_axis = time_axis/denominator

    num_pixels = dims[1]*dims[2]
    pixels = np.reshape(dataCube,(dims[0], num_pixels))
    pixels = pixels - np.mean(pixels, axis=0)

    slopeImage = np.dot(time_axis, pixels )
    slopeImage = np.reshape(slopeImage,(dims[1], dims[2]))

    #print('test slope image pixel at 110, 115: ', slopeImage[110,115])

# Slope image calculation
# calculated from equation in https://en.wikipedia.org/wiki/Simple_linear_regression
# mask saturated pixels
def calc_masked_slopeImage():
    print('in calc_masked_slopeImage()')
    global slopeImage, time_axis
    
    num_pixels = dims[1]*dims[2]
    pixels = np.reshape(dataCube,(dims[0], num_pixels))
    mask = (pixels == 65535)
    pixels = np.ma.masked_array(pixels, mask)
    num_masked = np.ma.count_masked(pixels)
    print('number of masked pixel elements: ', num_masked )
    
    if (num_masked > 0) :
  
        pixel_means = np.ma.mean(pixels, axis=0)
        pixels = pixels - pixel_means
        print('size of pixels: ', pixels.shape)


        # calculate denominator for masked elements
        #denominator = np.sum(time_axis**2)

        time_axis = np.ma.masked_array((np.repeat(time_axis.reshape(-1,1), num_pixels,1)),mask)
        time_means = np.ma.mean(time_axis, axis=0)
        time_axis = time_axis - time_means
        denominator =np.ma.sum(time_axis**2,axis=0)

        slopeImage = np.ma.sum(time_axis*pixels/denominator,axis=0)
        
    else : 
        print('no pixels needed to be masked')
        calc_raw_slopeImage()
    
#-----------------------------------
def write_level1b_output():
    global root_attrs, frames_attrs, outFile
    # In[7]: create output level1-b filename, from reference frame's acquisition date
    #filedate= hdu[0].header.get('DATE');
    
    #filedate = filedate.decode()
   # print('filedate: ', filedate)
   # filedate = os.path.splitext(filedate)[0]

   # year,month, day_hour = filedate.split('-')
   # day, hour = day_hour.split(' ')
   # hour, mins, seconds = hour.split(':')

    # get data, hour, minute image was taken
   # outputFile = outputPath+'/wfi_1b_'+year+month+day+hour+mins+'_00.h5'
    
    # write datasets
    if compression_flag :
        cds_dataset = outFile.create_dataset('cdsImage',  (dims[1],dims[2]), dtype='uint16',chunks=(128,128),compression='gzip', compression_opts=9, data=cdsImage)
        slope_dataset = outFile.create_dataset('slopeImage',(dims[1],dims[2]),dtype='float32',chunks=(128,128),compression='gzip', compression_opts=9, data=slopeImage)
    else :
        cds_dataset = outFile.create_dataset('cdsImage',  (dims[1],dims[2]), dtype='uint16',chunks=(128,128), data=cdsImage)
        slope_dataset= outFile.create_dataset('slopeImage',(dims[1],dims[2]),dtype='float32',chunks=(128,128), data=slopeImage)
   
    # copy input file root attributes to output file root attributes
    root_attrs = dict(root_attrs.items())
    outFile.attrs.update(root_attrs)
    
    
    # copy Frames dataset attributes to slope_dataset and cds_dataset attributes
    frames_attrs = dict(frames_attrs.items())
    frames_attrs['test_name'] = '-999'
    frames_attrs['test_phase'] = '-999'
    frames_attrs['test_script'] = '-999'
    frames_attrs['calibration_software'] = '-999'
    frames_attrs['calibration_version'] = np.uint16(0)
    frames_attrs['calibration_file'] = '-999'
    frames_attrs['filename_level1a'] = inputFile
    frames_attrs['filename'] = outputFile

    
    # copy frames attrs to the datasets
    cds_dataset.attrs.update(frames_attrs)
    slope_dataset.attrs.update(frames_attrs)
        
    #create Frames and Telementry links
    #g1 = outFile.create_group('/Frames')
    #g2 = outFile.create_group('/Telemetry')
    
    # check if Frames data/link exists

    print('/Frames: ')
    try: 
     frames = inFile["Frames"]
     print ('The Frames data exists in input file, creating an external link to it.\n')
     outFile["Frames"] = h5py.ExternalLink(inputFile,"/Frames")
    except: 
     print('The Frames data does NOT exist\n or some other error occurred.\n')
    
  
    print('/Telemetry: ')
    try: 
     telemetry = inFile["Telemetry"]
     print ('The Telemetry group exists in input file, creating an external link to it.\n')
     outFile["Telemetry"] = h5py.ExternalLink(inputFile,"/Telemetry")
    except: 
     print('The Telemetry group does NOT exist\n or could not create external link to it.\n')    
    
    #-------------------
    

    # currently provide the path to original .h5 Frames and Telemetry , update accordingly
  

    #close the file
    outFile.close()
    inFile.close() 
    
#-------------------------------------------------
# check command line parameters
try:
 inputFile = str(sys.argv[1])
except:
    print('Please run the code as follows: python3 level_1b_processing.py input_filename output_filename.h5 mask_flag compression_flag\n')
    sys.exit("Error: Enter Complete Path for the Input Filename")

if (inputFile.endswith(".h5")) :
    print('Input Level 1A .h5 filename: ', inputFile)
else:
    sys.exit("Extension for Level 1A input must be .h5")

# check being able to read from and write to input/output files
try:
 inFile = h5py.File(inputFile,'r')
except:
 print('Invalid input file path or name: '+inputFile)
 sys.exit('Please check the input file path and name, read permissions, and that the Level 1A input filename has .h5 extension')

try:
 outputFile = str(sys.argv[2])
except:
 print('Please run the code as follows: python3 level_1b_processing.py input_filename output_filename mask_flag compression_flag\n')
 sys.exit("Error: Enter complete path for Level 1b Output Filename ")
 
if (outputFile.endswith(".h5")) :
    print('Output file: ', outputFile)
else:
    sys.exit("Extension for Level 1B output must be .h5")

# create and write level-1b output file in the passed output directory
try:
 outFile = h5py.File(outputFile,'w')
except:
 print('Invalid output file path or name: '+outputFile)
 sys.exit('Please check the output file path and name, write permissions, and that the Level 1B output filename has .h5 extension\n')
    
try:
 mask_flag = int(sys.argv[3])
 print('mask_flag: ', mask_flag)
except:
 print('Please run the code as follows: python3 level_1b_processing.py input_filename.fits output_filename mask_flag compression_flag\n')   
 sys.exit("Error: Enter Mask_Flag (0 or 1) for Level 1b Output File ")
    
try:
 compression_flag = int(sys.argv[4])
 print('compression_flag: ', compression_flag)
except:
 print('Please run the code as follows: python3 level_1b_processing.py input_filename.fits output_filename mask_flag compression_flag\n')
 sys.exit("Error: Enter Compression_Flag (0 or 1) for Level 1b Output File ")

# Run the routines and time performance of the routine
#curr_time = time.process_time()

# readRSTFitsFile()
readRST_L1A_File()
calc_cdsImage()
if mask_flag :
    calc_masked_slopeImage()
else :
    calc_raw_slopeImage()
write_level1b_output()

print('End of Level 1b Processing.')
#elapsed_time = time.process_time() - curr_time
#print('Elapsed time for level-1b processing and writing the output: ',elapsed_time)


