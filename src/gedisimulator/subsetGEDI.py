
'''
Script to access simuated GEDI data
'''

##################################
import numpy as np
import h5py
from sys import exit
import matplotlib.pyplot as plt
if __name__ == '__main__':
  import argparse

###################################

class gediData(object):
  '''
  Simulated GEDI data handler
  '''

  def __init__(self,filename=None,minX=-100000000,maxX=100000000,minY=-1000000000,maxY=100000000,outName="teast.h5"):
    '''
    Class initialiser. Calls a function
    to read waveforms between bounds
    and writes to a new file
    '''

    self.subsetGEDI(filename,minX,maxX,minY,maxY,outName)


  ###########################################

  def subsetGEDI(self,filename,minX,maxX,minY,maxY,outNamen):
    '''
    Read real GEDI data from file
    '''
    # open file for reading
    f=h5py.File(filename,'r')
    self.beamList=['BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011']
    self.nWaves=0

    # open output file
    outFile=h5py.File(outNamen,'w')

    # set directory list
    self.setRealList()

    # loop over beams
    for b in self.beamList:
      if((b in list(f))==False): # does this exist?
        continue                 # if not, skip it
      elif(('geolocation' in list(f[b]))==False):  # no data in bea,
        continue

      print(b)

      # read the coords and determine output
      allLat=(np.array(f[b]['geolocation']['latitude_bin0'])+np.array(f[b]['geolocation']['latitude_lastbin']))/2.0
      allLon=(np.array(f[b]['geolocation']['longitude_bin0'])+np.array(f[b]['geolocation']['longitude_lastbin']))/2.0
      useInd=np.where((allLat>=minY)&(allLat<=maxY)&(allLon>=minX)&(allLon<=maxX))

      if(len(useInd[0])>0):
        useInd=useInd[0]
      else:      # none in here
        continue

      # create the beam group
      outFile.create_group(b)

      # loop over all arrays per shot
      for d in self.shotArrList:
        if((d=='rx_sample_start_index')|(d=='tx_sample_start_index')):  # skip these for noe
          continue
        # read array
        jimlad=np.array(f[b][d])[useInd]
        # write subset to a new file
        outFile[b].create_dataset(d,data=jimlad,compression='gzip')

      # for txwaveform and rxwaveform, we must read start/stop indices
      self.subsetWaves('rxwaveform','rx_sample_start_index','rx_sample_count',useInd,outFile[b],f[b])
      self.subsetWaves('txwaveform','tx_sample_start_index','tx_sample_count',useInd,outFile[b],f[b])

      # geolocation data
      g='geolocation'
      outFile[b].create_group(g)
      for d in self.geoArrList:
        if(d!='surface_type'):
          jimlad=np.array(f[b][g][d])[useInd]
          outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')
        else:
          jimlad=np.array(f[b][g][d])[:,useInd]
          outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')

      # ancillary data
      g='ancillary'
      outFile[b].create_group(g)
      for d in self.ancArrList:
        jimlad=np.array(f[b][g][d])
        outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')

      # geophys_corr
      g='geophys_corr'
      outFile[b].create_group(g)
      for d in self.corArrList:
        jimlad=np.array(f[b][g][d])[useInd]
        outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')

    f.close()
    outFile.close()
    print("Written to",outNamen)
    return


  ###########################################

  def subsetWaves(self,waveName,indName,lenName,useInd,outFile,f):
    '''Subset and write the RX or TX waveforms'''

    # read indices
    startInds=np.array(f[indName])[useInd]
    lenInds=np.array(f[lenName])[useInd]
    totBins=np.sum(lenInds)
    waveform=np.empty(totBins,dtype=np.float32)
    newInds=np.zeros(len(useInd),dtype=np.uint64)

    # read raw data and repack
    jimlad=np.array(f[waveName])
    lastInd=0
    for i in range(0,startInds.shape[0]):
      newInds[i]=lastInd
      waveform[lastInd:lastInd+lenInds[i]]=jimlad[startInds[i]:startInds[i]+lenInds[i]]
      lastInd=lastInd+lenInds[i]

    # write data
    outFile.create_dataset(waveName,data=waveform,compression='gzip')
    outFile.create_dataset(indName,data=newInds,compression='gzip')

    return


  ###########################################

  def setRealList(self):
    '''Set list of all data within a real GEDI file'''

    # arrays with one element per shot
    self.shotArrList=['all_samples_sum', 'beam', 'channel', 'delta_time', 'master_frac', 'master_int',\
                      'noise_mean_corrected', 'noise_stddev_corrected', 'nsemean_even', 'nsemean_odd',\
                      'rx_energy', 'rx_offset', 'rx_open', 'rx_sample_count',\
                      'selection_stretchers_x', 'selection_stretchers_y', 'shot_number', 'stale_return_flag',\
                      'th_left_used', 'tx_egamplitude', 'tx_egamplitude_error', 'tx_egbias', 'tx_egbias_error',\
                      'tx_egflag', 'tx_eggamma', 'tx_eggamma_error', 'tx_egsigma', 'tx_egsigma_error', 'tx_gloc',\
                      'tx_gloc_error', 'tx_pulseflag', 'tx_sample_count',]

    self.geoArrList=['altitude_instrument', 'altitude_instrument_error', 'bounce_time_offset_bin0', 'bounce_time_offset_bin0_error',\
                     'bounce_time_offset_lastbin', 'bounce_time_offset_lastbin_error', 'degrade', 'delta_time',\
                     'digital_elevation_model', 'elevation_bin0', 'elevation_bin0_error', 'elevation_lastbin',\
                     'elevation_lastbin_error', 'latitude_bin0', 'latitude_bin0_error', 'latitude_instrument',\
                     'latitude_instrument_error', 'latitude_lastbin', 'latitude_lastbin_error', 'local_beam_azimuth',\
                     'local_beam_azimuth_error', 'local_beam_elevation', 'local_beam_elevation_error', 'longitude_bin0',\
                     'longitude_bin0_error', 'longitude_instrument', 'longitude_instrument_error', 'longitude_lastbin',\
                     'longitude_lastbin_error', 'mean_sea_surface', 'neutat_delay_derivative_bin0', 'neutat_delay_derivative_lastbin',\
                     'neutat_delay_total_bin0', 'neutat_delay_total_lastbin', 'range_bias_correction', 'shot_number',\
                     'solar_azimuth', 'solar_elevation', 'surface_type']

    self.ancArrList=['master_time_epoch', 'mean_samples', 'smoothing_width']

    self.corArrList=['delta_time', 'dynamic_atmosphere_correction', 'geoid', 'tide_earth', 'tide_load', 'tide_ocean',\
                     'tide_ocean_pole', 'tide_pole']

    return

  ###########################################

  def findBounds(self,meanN,stdev,i):
    '''Find the signal start and end'''
    thresh=3.5*stdev+meanN
    minWidth=3
    binList=np.where(self.wave[i]>thresh)
    buff=15

    topBin=0
    for j in range(0,len(binList[0])):
      if (binList[0][j]==(binList[0][j-1]+1))&(binList[0][j]==(binList[0][j-2]+2)):
        topBin=binList[0][j]
        break

    botBin=binList[len(binList)-1]
    for j in range(len(binList[0])-1,0,-1):
      if (binList[0][j]==(binList[0][j-1]+1))&(binList[0][j]==(binList[0][j-2]+2)):
        botBin=binList[0][j]
        break

    return(self.z[botBin]-buff,self.z[topBin]+buff)

  ###########################################

  def writeCoords(self):
    for i in range(0,len(self.lon)):
      print(self.lon[i],self.lat[i])


# end of gediData class
###########################################


###########################################
# read the command line

if __name__ == '__main__':
  def gediCommands():
    '''
    Read commandline arguments
    '''
    p = argparse.ArgumentParser(description=("Writes out properties of GEDI waveform files"))
    p.add_argument("--input",dest="inName",type=str,help=("Input GEDI HDF5 filename"))
    p.add_argument("--bounds", dest ="bounds", type=float,nargs=4,default=[-100000000,-100000000,100000000000,10000000000], help=("Bounds to plot between. minX minY maxX maxY"))
    p.add_argument("--output",dest="output",type=str,default='teast.h5',help=("Output filename"))
    cmdargs = p.parse_args()
    return cmdargs


###########################################
# the main block

if __name__ == '__main__':
  # read the command line
  cmdargs=gediCommands()
  inName=cmdargs.inName
  bounds=cmdargs.bounds
  output=cmdargs.output

  # read data
  gedi=gediData(filename=inName,minX=bounds[0],maxX=bounds[2],minY=bounds[1],maxY=bounds[3],outName=output)

