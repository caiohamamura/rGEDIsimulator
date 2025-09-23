#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "inttypes.h"
#include "tools.h"
#include "functionWrappers.h"
#include "hdf5.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "libLasRead.h"
#include "libDEMhandle.h"
#include "libLidVoxel.h"
#include "libOctree.h"
#include "gediIO.h"
#include "gediNoise.h"
#ifndef WITHOUT_GDAL
  #include "ogr_srs_api.h"
#endif
#include "gsl/gsl_fft_complex.h"


/*tolerances*/
#define TOL 0.00001



/*####################################################*/
/*tidy data structure*/

dataStruct **tidyAsciiStruct(dataStruct **data,int nFiles)
{
  int i=0;

  if(data){
    for(i=0;i<nFiles;i++){
      TTIDY((void **)data[i]->wave,data[i]->nWaveTypes);
      TTIDY((void **)data[i]->ground,data[i]->nWaveTypes);
      TIDY(data[i]->noised);
      TIDY(data[i]->totE);
      TIDY(data[i]->z);
      TIDY(data[i]);
    }
    TIDY(data);
  }

  return(data);
}/*tidyAsciiStruct*/


/*####################################################*/
/*read ASCII data*/

dataStruct *readASCIIdata(char *namen,gediIOstruct *gediIO)
{
  int i=0,ind=0,*numb=NULL;
  dataStruct *data=NULL;
  char line[1000],temp1[100],temp2[100];
  char temp3[100],temp4[100],temp5[100];
  char temp6[100],temp7[100],temp8[100];
  char temp9[100],temp10[100];
  FILE *ipoo=NULL;

  /*open input*/
  if((ipoo=fopen(namen,"r"))==NULL){
    errorf("Error opening input file %s\n",namen);
    return(NULL);
  }

  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    errorf("error control allocation.\n");
    return(NULL);
  }
  data->pSigma=-1.0;    /*nonesense pulse length*/
  data->fSigma=-1.0;    /*nonesense footprint width*/
  data->usable=1;
  data->zen=0.0;        /*straight down*/
  data->nWaveTypes=(int)(gediIO->useInt+gediIO->useCount+gediIO->useFrac);
  if(data->nWaveTypes==1)data->useType=0;

  /*count number of wavebins*/
  data->nBins=0;
  while(fgets(&(line[0]),1000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1))data->nBins++;
  }


  /*is there usable data?*/
  if(data->nBins==0){
    data->usable=0;
  }else{
    ASSIGN_CHECKNULL_RETNULL(data->wave,fFalloc(data->nWaveTypes,"waves",0));
    if(gediIO->ground) {ASSIGN_CHECKNULL_RETNULL(data->ground,fFalloc(data->nWaveTypes,"ground",0));}
    ASSIGN_CHECKNULL_RETNULL(data->z,dalloc(data->nBins,"z",0));
    for(i=0;i<data->nWaveTypes;i++){
    ASSIGN_CHECKNULL_RETNULL(data->wave[i],falloc((uint64_t)data->nBins,"waveform",0));
      if(gediIO->ground) {ASSIGN_CHECKNULL_RETNULL(data->ground[i],falloc((uint64_t)data->nBins,"ground",0));}
    }
    data->useID=0;
    data->demGround=0;

    /*rewind to start of file*/
    if(fseek(ipoo,(long)0,SEEK_SET)){
      errorf("fseek error\n");
      return(NULL);
    }

    /*read data*/
    ASSIGN_CHECKNULL_RETNULL(numb,ialloc(data->nWaveTypes,"number",0));
    for(i=0;i<data->nWaveTypes;i++)numb[i]=0;
    while(fgets(line,800,ipoo)!=NULL){
      if(strncasecmp(line,"#",1)){
        if(gediIO->useInt){  /*read intensity*/
          ind=0;
          if(gediIO->ground==0){   /*don't read ground*/
            if(sscanf(line,"%s %s",temp1,temp2)==2){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp2);
              numb[ind]++;
            }
          }else{                   /*read ground*/
            if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp2);
              data->ground[ind][numb[ind]]=atof(temp4);
              numb[ind]++;
            }
          }/*ground switch*/
        }
        if(gediIO->useFrac){ /*read fraction*/
          ind=(int)(gediIO->useInt+gediIO->useCount);
          if(gediIO->ground==0){   /*don't read ground*/
           if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp3);
              numb[ind]++;
            }
          }else{                   /*read ground*/
            if(sscanf(line,"%s %s %s %s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10)==10){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp8);
              data->ground[ind][numb[ind]]=atof(temp10);
              numb[ind]++;
            }
          }/*ground switch*/
        }
        if(gediIO->useCount){ /*read count*/
          ind=(int)gediIO->useInt;
          if(gediIO->ground==0){   /*don't read ground*/
            if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp3);
              numb[ind]++;
            }
          }else{                   /*read ground*/
            if(sscanf(line,"%s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7)==7){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp5);
              data->ground[ind][numb[ind]]=atof(temp7);
              numb[ind]++;
            }
          }/*ground switch*/
        }/*intensity or count switch*/
      }else{  /*read the header*/
        if(sscanf(line,"%s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7)==7){
          if(!strncasecmp(temp2,"fSigma",6)){
            data->fSigma=atof(temp3);
            data->pSigma=atof(temp5);
          }
        }
        if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
          if(!strncasecmp(temp2,"waveID",6)){
            data->useID=1;
            strcpy(&(data->waveID[0]),temp3);
          }else if(!strncasecmp(temp2,"meanScanAng",11)){
            data->zen=atof(temp3);
          }
        }
        if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
          if(!strncasecmp(temp2,"coord",5)){
            data->lon=atof(temp3);
            data->lat=atof(temp4);
          }else if(!strncasecmp(temp2,"ground",6)){
            if(!gediIO->dontTrustGround){
              data->gElev=atof(temp3);
              data->slope=atof(temp4);
              data->demGround=1;
            }else data->demGround=0;
          }
        }
        if(sscanf(line,"%s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6)==6){
          if(!strncasecmp(temp2,"density",7)){
            data->pointDense=atof(temp4);
            data->beamDense=atof(temp6);
          }else if(!strncasecmp(temp2,"lvis",4)){
            data->res=atof(temp4);
            data->zen=atof(temp6);
          }
        }
      }
    }/*line loop*/
    TIDY(numb);

    /*add up energy*/
    ASSIGN_CHECKNULL_RETNULL(data->totE,falloc((uint64_t)data->nWaveTypes,"",0));
    for(ind=0;ind<data->nWaveTypes;ind++){
      data->totE[ind]=0.0;
      for(i=0;i<data->nBins;i++)data->totE[ind]+=data->wave[ind][i];
    }

    gediIO->res=gediIO->den->res=gediIO->gFit->res=fabs(data->z[data->nBins-1]-data->z[0])/(float)(data->nBins-1);
    if(data->res<=0.0)data->res=gediIO->res;
    if(gediIO->den->res<TOL)data->usable=0;
    if(data->totE[data->useType]<=0.0)data->usable=0;
    if(gediIO->ground==0){   /*set to blank*/
      data->cov=-1.0;
      data->gLap=-1.0;
      data->gMinimum=-1.0;
      data->gInfl=-1.0;
    }
    if(data->fSigma<0.0)data->fSigma=gediIO->fSigma;
    if(data->pSigma<0.0)data->pSigma=gediIO->pSigma;
  }/*check there is dara*/

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  /*set up number of messages*/
  if((gediIO->nFiles>gediIO->nMessages)&&(gediIO->nMessages>0))gediIO->nMessages=(int)(gediIO->nFiles/gediIO->nMessages);
  else                                gediIO->nMessages=1;

  return(data);
}/*readASCIIdata*/


/*####################################################*/
/*turn data into HDF5 structure*/

gediHDF *arrangeGEDIhdf(dataStruct **data,gediIOstruct *gediIO)
{
  int i=0;
  gediHDF *hdfData=NULL;
  int trimDataLength(dataStruct **,gediHDF *,gediIOstruct *);

  /*allocate space for all*/
  if(!(hdfData=(gediHDF *)calloc(1,sizeof(gediHDF)))){
    errorf("error control allocation.\n");
    return(NULL);
  }

  /*count nuber of usable waves*/
  hdfData->nWaves=0;
  for(i=0;i<gediIO->nFiles;i++)if(data[i]->usable)hdfData->nWaves++;

  /*per scan*/
  hdfData->pSigma=data[0]->pSigma;
  hdfData->fSigma=data[0]->fSigma;
  hdfData->nPbins=0;
  hdfData->pulse=NULL;
  /*per beam*/
  ASSIGN_CHECKNULL_RETNULL(hdfData->z0,falloc((uint64_t)hdfData->nWaves,"top elevations",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->zN,falloc((uint64_t)hdfData->nWaves,"bottom elevations",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->lon,dalloc(hdfData->nWaves,"lon",0));     /*longitudes*/
  ASSIGN_CHECKNULL_RETNULL(hdfData->lat,dalloc(hdfData->nWaves,"lat",0));    /*latitudes*/
  ASSIGN_CHECKNULL_RETNULL(hdfData->slope,falloc((uint64_t)hdfData->nWaves,"slope",0));    /*ground slope*/
  ASSIGN_CHECKNULL_RETNULL(hdfData->gElev,falloc((uint64_t)hdfData->nWaves,"ground elevation, CofG",0));    /*ground elevation, CofG*/
  ASSIGN_CHECKNULL_RETNULL(hdfData->demElev,falloc((uint64_t)hdfData->nWaves,"ground elevation, DEM",0));  /*ground elevation, DEM*/
  ASSIGN_CHECKNULL_RETNULL(hdfData->beamDense,falloc((uint64_t)hdfData->nWaves,"beamDense",0));/*beam density*/
  ASSIGN_CHECKNULL_RETNULL(hdfData->pointDense,falloc((uint64_t)hdfData->nWaves,"pointDense",0));/*point density*/
  ASSIGN_CHECKNULL_RETNULL(hdfData->zen,falloc((uint64_t)hdfData->nWaves,"zen",0));      /*scan angles, or mean angles*/


  /*trim and copy data*/
  ISINTRETNULL(trimDataLength(data,hdfData,gediIO));

  return(hdfData);
}/*arrangeGEDIhdf*/


/*####################################################*/
/*trim all arrays to be the same length*/

int trimDataLength(dataStruct **data,gediHDF *hdfData,gediIOstruct *gediIO)
{
  int i=0,j=0,maxBins=0,maxID=0;
  int ind=0,numb=0;;
  uint64_t place=0;
  int *start=NULL,*end=NULL;
  float tot=0,thresh=0,res=0;
  float buffer=0,cumul=0;

  /*buffer length in metres*/
  if(gediIO->pcl==0)buffer=30.0;
  else              buffer=0.0;

  /*allocate usable bins*/
  ASSIGN_CHECKNULL_RETINT(start,ialloc(hdfData->nWaves,"wave starts",0));
  ASSIGN_CHECKNULL_RETINT(end,ialloc(hdfData->nWaves,"end starts",0));


  /*determine start and end of all waves*/
  maxBins=maxID=-1;
  for(i=0;i<hdfData->nWaves;i++){
    if(!data[i]->usable)continue;
    for(ind=0;ind<data[i]->nWaveTypes;ind++){
      /*total energy for a threshod*/
      tot=0.0;
      for(j=0;j<data[i]->nBins;j++)tot+=data[i]->wave[ind][j];
      thresh=0.0005*tot;

      /*determine used bounds*/
      cumul=0.0;
      for(j=0;j<data[i]->nBins;j++){
        cumul+=data[i]->wave[ind][j];
        if(cumul>=thresh){
          start[i]=j;
          break;
        }
      }
      cumul=0.0;
      for(j=data[i]->nBins-1;j>=0;j--){
        cumul+=data[i]->wave[ind][j];
        if(cumul>=thresh){
          end[i]=j;
          break;
        }
      }
    }

    /*add a buffer for later smoothing*/
    res=fabs(data[i]->z[0]-data[i]->z[data[i]->nBins-1])/(float)(data[i]->nBins-1);
    start[i]-=(int)(buffer/res);
    end[i]+=(int)(buffer/res);
    if(start[i]<0)start[i]=0;
    if(end[i]>=data[i]->nBins)end[i]=data[i]->nBins-1;

    /*determine max*/
    if((end[i]-start[i])>maxBins)maxBins=end[i]-start[i];
    if(data[i]->useID){
      if(((int)strlen(data[i]->waveID)+1)>maxID)maxID=(int)strlen(data[i]->waveID)+1;
    }
  }

  /*if we are doing PCL, do not zero pad and just save the pulse*/
  if(gediIO->pcl+gediIO->pclPhoton)maxBins=(int)((float)hdfData->nPbins*hdfData->pRes/gediIO->res*4.0);

  ASSIGN_CHECKNULL_RETINT(hdfData->nBins,ialloc(1,"bins",0));
  hdfData->nBins[0]=maxBins;

  if(maxID>0)hdfData->idLength=maxID;
  else       hdfData->idLength=7;

  /*allocate funny arrays needed by HDF5 library*/
  hdfData->nTypeWaves=data[0]->nWaveTypes;
  ASSIGN_CHECKNULL_RETINT(hdfData->wave,fFalloc(hdfData->nTypeWaves,"waveforms",0));
  ASSIGN_CHECKNULL_RETINT(hdfData->ground,fFalloc(hdfData->nTypeWaves,"ground waves",0));
  for(ind=0;ind<hdfData->nTypeWaves;ind++){
    ASSIGN_CHECKNULL_RETINT(hdfData->wave[ind],falloc((uint64_t)hdfData->nWaves*(uint64_t)hdfData->nBins[0],"waveforms",0));
    ASSIGN_CHECKNULL_RETINT(hdfData->ground[ind],falloc((uint64_t)hdfData->nWaves*(uint64_t)hdfData->nBins[0],"ground waves",0));
  }
  ASSIGN_CHECKNULL_RETINT(hdfData->waveID,challoc(hdfData->nWaves*hdfData->idLength,"wave IDs",0));

  /*copy arrays*/
  numb=0;
  for(i=0;i<hdfData->nWaves;i++){
    if(!data[i]->usable)continue;
    /*range and resolution*/
    res=fabs(data[i]->z[0]-data[i]->z[data[i]->nBins-1])/(float)(data[i]->nBins-1);
    hdfData->z0[numb]=data[i]->z[start[i]];
    hdfData->zN[numb]=hdfData->z0[numb]-(float)hdfData->nBins[0]*res;
    /*copy data*/
    for(ind=0;ind<hdfData->nTypeWaves;ind++){
      for(j=start[i];j<end[i];j++){
        place=(uint64_t)numb*(uint64_t)hdfData->nBins[0]+(uint64_t)j-(uint64_t)start[i];
        hdfData->wave[ind][place]=data[i]->wave[ind][j];
        hdfData->ground[ind][place]=data[i]->ground[ind][j];
      }
      /*pad end if not long enough*/
      for(j=end[i];j<(maxBins+start[i]);j++){
        place=(uint64_t)numb*(uint64_t)hdfData->nBins[0]+(uint64_t)j-(uint64_t)start[i];
        hdfData->wave[ind][place]=0.0;
        hdfData->ground[ind][place]=0.0;
      }
    }

    hdfData->lon[numb]=data[i]->lon;
    hdfData->lat[numb]=data[i]->lat;
    hdfData->slope[numb]=data[i]->slope;
    hdfData->gElev[numb]=data[i]->gElev;
    hdfData->demElev[numb]=data[i]->demGround;
    if(maxID>0)strcpy(&(hdfData->waveID[numb*hdfData->idLength]),data[i]->waveID);
    else       snprintf(&(hdfData->waveID[numb*hdfData->idLength]), hdfData->idLength,"%d",i);
    hdfData->beamDense[numb]=data[i]->beamDense;
    hdfData->pointDense[numb]=data[i]->pointDense;
    hdfData->zen[numb]=data[i]->zen;
    numb++;
  }/*wave loop*/

  TIDY(start);
  TIDY(end);
  return(0);
}/*trimDataLength*/


/*####################################################*/
/*write data to HDF5 in L1B format*/

int writeGEDIl1b(gediHDF *hdfData,char *namen,gediIOstruct *gediIO)
{
  #ifndef WITHOUT_GDAL
  int i=0;
  uint8_t *tempUint8=NULL;
  uint8_t *padUint8zeros(int);
  uint8_t *padUint8ones(int);
  uint16_t *tempUint16=NULL;
  uint16_t *padUint16zeros(int);
  uint16_t *setRxSampleCount(int *,int);
  uint16_t *setSelectStretchL1B(int);
  uint16_t *setThUsedL1B(int);
  uint32_t *tempUint32=NULL;
  uint32_t *padUint32zeros(int);
  uint32_t *setAllSamplesSumL1B(gediHDF *);
  uint64_t *tempUint64=NULL;
  uint64_t *setRXstarts(int,int *);
  uint64_t *setShotNumber(int);
  int8_t *tempInt8=NULL;
  int8_t *setSurfaceTypeL1B(int,int);
  int32_t *padInt32ones(int);
  int64_t *tempInt64=NULL;
  float *tempFloat=NULL;
  float *setTXegAmpL1B(int,float);
  double *tempDouble=NULL;
  double *setAltitude(int);
  double *setBounceOffset(int,int *,float);
  double *setDeltaTime(int);
  double *setDEMl1b(gediHDF *,int);
  double *setElevBinL1B(int,float *);
  double *setCovL1B(gediHDF *,int);
  double *setL1Bcoords(int,gediHDF *);
  double *setHalfPiL1B(int);
  double *setDelayDerivL1B(int);
  double *setRXenergyL1B(gediHDF *);
  hid_t file,group_id,sgID;         /* Handles */
  herr_t status;
  TXstruct tx;          /*to hold pulse information for TX*/
  int rearrangePulsetoTX(gediIOstruct *,gediHDF *,TXstruct *);

  /*open new file*/
  file=H5Fcreate(namen,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  /*create the beam group*/
  group_id=H5Gcreate2(file,"BEAM0000", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /*write the data*/
  ASSIGN_CHECKNULL_RETINT(tempUint32,setAllSamplesSumL1B(hdfData));
  ISINTRETINT(writeComp1dUint32HDF5(group_id,"all_samples_sum",tempUint32,hdfData->nWaves));
  TIDY(tempUint32);
  ASSIGN_CHECKNULL_RETINT(tempUint16, padUint16zeros(hdfData->nWaves));   /*padded zeroes for fake beams*/
  ISINTRETINT(writeComp1dUint16HDF5(group_id,"beam",tempUint16,hdfData->nWaves));
  TIDY(tempUint16);
  ASSIGN_CHECKNULL_RETINT(tempUint8, padUint8zeros(hdfData->nWaves));   /*padded zeroes for fake beams*/
  ISINTRETINT(writeComp1dUint8HDF5(group_id,"channel",tempUint8,hdfData->nWaves));
  TIDY(tempUint8);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setDeltaTime(hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(group_id,"delta_time",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempFloat, falloc(hdfData->nWaves,"tempFloat",0));
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"master_frac",tempFloat,hdfData->nWaves));
  TIDY(tempFloat);
  ASSIGN_CHECKNULL_RETINT(tempUint32, padUint32zeros(hdfData->nWaves));   /*padded zeroes for fake beams*/
  ISINTRETINT(writeComp1dUint32HDF5(group_id,"master_int",tempUint32,hdfData->nWaves));
  TIDY(tempUint32);
  ASSIGN_CHECKNULL_RETINT(tempFloat, falloc(hdfData->nWaves,"tempFloat",0));    /*needs updating if noise added*/
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"mean",tempFloat,hdfData->nWaves));
  TIDY(tempFloat);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"noise_mean_corrected",0));
  ISINTRETINT(writeComp1dDoubleHDF5(group_id,"noise_mean_corrected",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"noise_stddev_corrected",0));
  ISINTRETINT(writeComp1dDoubleHDF5(group_id,"noise_stddev_corrected",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"nsemean_even",0));
  ISINTRETINT(writeComp1dDoubleHDF5(group_id,"nsemean_even",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"nsemean_odd",0));
  ISINTRETINT(writeComp1dDoubleHDF5(group_id,"nsemean_odd",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setRXenergyL1B(hdfData));
  ISINTRETINT(writeComp1dDoubleHDF5(group_id,"rx_energy",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempUint16, padUint16zeros(hdfData->nWaves));
  ISINTRETINT(writeComp1dUint16HDF5(group_id,"rx_offset",tempUint16,hdfData->nWaves));
  TIDY(tempUint16);
  ASSIGN_CHECKNULL_RETINT(tempUint32, padUint32zeros(hdfData->nWaves));
  ISINTRETINT(writeComp1dUint32HDF5(group_id,"rx_open",tempUint32,hdfData->nWaves));
  TIDY(tempUint32);
  ASSIGN_CHECKNULL_RETINT(tempUint16, setRxSampleCount(hdfData->nBins,hdfData->nWaves));
  ISINTRETINT(writeComp1dUint16HDF5(group_id,"rx_sample_count",tempUint16,hdfData->nWaves));
  TIDY(tempUint16);
  ASSIGN_CHECKNULL_RETINT(tempUint64,setRXstarts(hdfData->nWaves,hdfData->nBins));
  ISINTRETINT(writeComp1dUint64HDF5(group_id,"rx_sample_start_index",tempUint64,hdfData->nWaves));
  TIDY(tempUint64);
  if(gediIO->useCount){
    ISINTRETINT(writeComp1dFloatHDF5(group_id,"rxwaveform",hdfData->wave[(int)gediIO->useInt],hdfData->nWaves*hdfData->nBins[0]));
    if(hdfData->ground) { ISINTRETINT(writeComp1dFloatHDF5(group_id,"grxwaveform",hdfData->ground[(int)gediIO->useInt],hdfData->nWaves*hdfData->nBins[0])); }
  }else{
    errorf("Issues with HDF5 format and not using the count method\n");
    return -1;
  }
  ASSIGN_CHECKNULL_RETINT(tempUint16,setSelectStretchL1B(hdfData->nWaves));
  ISINTRETINT(writeComp1dUint16HDF5(group_id,"selection_stretchers_x",tempUint16,hdfData->nWaves));
  ISINTRETINT(writeComp1dUint16HDF5(group_id,"selection_stretchers_y",tempUint16,hdfData->nWaves));
  TIDY(tempUint16);
  ASSIGN_CHECKNULL_RETINT(tempUint64,setShotNumber(hdfData->nWaves));
  ISINTRETINT(writeComp1dUint64HDF5(group_id,"shot_number",tempUint64,hdfData->nWaves));
  TIDY(tempUint64);
  ASSIGN_CHECKNULL_RETINT(tempUint8, padUint8zeros(hdfData->nWaves));
  ISINTRETINT(writeComp1dUint8HDF5(group_id,"stale_return_flag",tempUint8,hdfData->nWaves));
  TIDY(tempUint8);
  ASSIGN_CHECKNULL_RETINT(tempUint16,setThUsedL1B(hdfData->nWaves));
  ISINTRETINT(writeComp1dUint16HDF5(group_id,"th_left_used",tempUint16,hdfData->nWaves));
  TIDY(tempUint16);
  ASSIGN_CHECKNULL_RETINT(tempFloat, falloc(hdfData->nWaves,"tempFloat",0));   /*needs updating if noise added*/
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"tx_egamplitude_error",tempFloat,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"tx_egbias",tempFloat,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"tx_egbias_error",tempFloat,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"tx_eggamma",tempFloat,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"tx_eggamma_error",tempFloat,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"tx_egsigma",tempFloat,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"tx_egsigma_error",tempFloat,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"tx_gloc",tempFloat,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"tx_gloc_error",tempFloat,hdfData->nWaves));
  TIDY(tempFloat);
  ASSIGN_CHECKNULL_RETINT(tempUint8, padUint8zeros(hdfData->nWaves));
  ISINTRETINT(writeComp1dUint8HDF5(group_id,"tx_egflag",tempUint8,hdfData->nWaves));
  TIDY(tempUint8);
  ASSIGN_CHECKNULL_RETINT(tempUint8, padUint8ones(hdfData->nWaves));
  ISINTRETINT(writeComp1dUint8HDF5(group_id,"tx_pulseflag",tempUint8,hdfData->nWaves));
  TIDY(tempUint8);
  ASSIGN_CHECKNULL_RETINT(tempFloat, falloc(hdfData->nWaves,"tempFloat",0));   /*needs updating if noise added*/
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"stddev",tempFloat,hdfData->nWaves));
  TIDY(tempFloat);

  /*rearrange the pulse in to the TX format*/
  ISINTRETINT(rearrangePulsetoTX(gediIO,hdfData,&tx));
  ISINTRETINT(writeComp1dUint16HDF5(group_id,"tx_sample_count",tx.txCount,hdfData->nWaves));
  ISINTRETINT(writeComp1dUint64HDF5(group_id,"tx_sample_start_index",tx.txStart,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"txwaveform",tx.txwave,hdfData->nWaves*(int)tx.nBins));
  ASSIGN_CHECKNULL_RETINT(tempFloat, setTXegAmpL1B(hdfData->nWaves,tx.maxAmp));
  ISINTRETINT(writeComp1dFloatHDF5(group_id,"tx_egamplitude",tempFloat,hdfData->nWaves));
  TIDY(tempFloat);


  /*tidy up*/
  TIDY(tx.txCount);
  TIDY(tx.txStart);
  TIDY(tx.txwave);

  /*add the sub-groups*/
  /*ancillary subgroup*/
  sgID=H5Gcreate2(group_id,"ancillary",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(1,"temp double",0));
  tempDouble[0]=(double)rand();
  ISINTRETINT(write1dDoubleHDF5(sgID,"master_time_epoch",tempDouble,1));
  TIDY(tempDouble);
  if(!(tempInt64=(int64_t *)calloc(1,sizeof(int64_t)))){
    errorf("error in tempInt64 allocation.\n");
    return -1;
  } 
  tempInt64[0]=100;
  ISINTRETINT(write1dInt64HDF5(sgID,"mean_samples",tempInt64,1));
  TIDY(tempInt64);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(1,"temp double",0));
  tempDouble[0]=0.0;
  ISINTRETINT(write1dDoubleHDF5(sgID,"smoothing_width",tempDouble,1));
  TIDY(tempDouble);
  status=H5Gclose(sgID);
 if(status<0){
    errorf("Error closing HDF5 group\n");
    return -1;
  }

  /*geolocation subgroup*/
  sgID=H5Gcreate2(group_id,"geolocation",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setAltitude(hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"altitude_instrument",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"altitude_instrument_error",0));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"altitude_instrument_error",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setBounceOffset(hdfData->nWaves,NULL,gediIO->res));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"bounce_time_offset_bin0",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"bounce_time_offset_bin0_error",0));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"bounce_time_offset_bin0_error",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setBounceOffset(hdfData->nWaves,hdfData->nBins,gediIO->res));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"bounce_time_offset_lastbin",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"bounce_time_offset_lastbin_error",0));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"bounce_time_offset_lastbin_error",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempUint8,padUint8zeros(hdfData->nWaves));
  ISINTRETINT(writeComp1dUint8HDF5(sgID,"degrade",tempUint8,hdfData->nWaves));
  TIDY(tempUint8);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setDeltaTime(hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"delta_time",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setDEMl1b(hdfData,(int)gediIO->useInt));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"digital_elevation_model",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setElevBinL1B(hdfData->nWaves,hdfData->z0));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"elevation_bin0",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"elevation_bin0_error",0));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"elevation_bin0_error",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setElevBinL1B(hdfData->nWaves,hdfData->zN));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"elevation_lastbin",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"elevation_lastbin_error",0));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"elevation_lastbin_error",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setCovL1B(hdfData,(int)gediIO->useInt));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"landsat_treecover",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"modis_treecover",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  /*reproject bounds to 4326 for GEDI*/
  ASSIGN_CHECKNULL_RETINT(tempDouble, setL1Bcoords(gediIO->aEPSG,hdfData));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"longitude_bin0",&(tempDouble[0]),hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"longitude_lastbin",&(tempDouble[0]),hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"longitude_instrument",&(tempDouble[0]),hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"latitude_bin0",&(tempDouble[hdfData->nWaves]),hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"latitude_lastbin",&(tempDouble[hdfData->nWaves]),hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"latitude_instrument",&(tempDouble[hdfData->nWaves]),hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"latitude_bin0_error",0));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"latitude_bin0_error",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"latitude_instrument_error",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"latitude_lastbin_error",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"local_beam_azimuth",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"local_beam_azimuth_error",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"local_beam_elevation_error",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"longitude_bin0_error",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"longitude_instrument_error",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"longitude_lastbin_error",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"mean_sea_surface",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"modis_nonvegetated",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"modis_nonvegetated_sd",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"modis_treecover_sd",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setHalfPiL1B(hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"local_beam_elevation",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setDelayDerivL1B(hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"neutat_delay_derivative_bin0",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"neutat_delay_derivative_lastbin",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"neutat_delay_total_bin0",0));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"neutat_delay_total_bin0",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"neutat_delay_total_lastbin",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"range_bias_correction",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"solar_azimuth",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"solar_elevation",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempUint64,setShotNumber(hdfData->nWaves));
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"shot_number",0));
  for(i=0;i<hdfData->nWaves;i++)tempDouble[i]=(double)tempUint64[i];
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"shot_number",tempDouble,hdfData->nWaves));
  TIDY(tempUint64);
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempInt8,setSurfaceTypeL1B(hdfData->nWaves,5));
  ISINTRETINT(writeComp2dInt8HDF5(sgID,"surface_type",tempInt8,5,hdfData->nWaves));
  TIDY(tempInt8);
  status=H5Gclose(sgID);
 if(status<0){
    errorf("Error closing HDF5 group\n");
    return -1;
  }

  /*geophys_corr subgroup*/
  sgID=H5Gcreate2(group_id,"geophys_corr",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  ASSIGN_CHECKNULL_RETINT(tempDouble, setDeltaTime(hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"delta_time",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  ASSIGN_CHECKNULL_RETINT(tempDouble, dalloc(hdfData->nWaves,"dynamic_atmosphere_correction",0));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"dynamic_atmosphere_correction",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"geoid",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"tide_earth",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"tide_load",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"tide_ocean",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"tide_ocean_pole",tempDouble,hdfData->nWaves));
  ISINTRETINT(writeComp1dDoubleHDF5(sgID,"tide_pole",tempDouble,hdfData->nWaves));
  TIDY(tempDouble);
  status=H5Gclose(sgID);
  if(status<0){
    errorf("Error closing HDF5 group\n");
    return -1;
  }


  /*close the beam group*/
  status=H5Gclose(group_id);
 if(status<0){
    errorf("Error closing HDF5 group\n");
    return -1;
  }

  /*close file*/
  if(H5Fclose(file)){
    errorf("Issue closing file\n");
    return -1;
  }

  msgf("Waveforms written to %s\n",namen);
  return 0;
#else
  return 0;
#endif //WITHOUT_GDAL
}/*writeGEDIl1b*/


/*####################################################*/
/*populate array with  TX max amplitude*/

float *setTXegAmpL1B(int nWaves,float maxAmp)
{
  int i=0;
  float *tempFloat=NULL;

  ASSIGN_CHECKNULL_RETNULL(tempFloat,falloc(nWaves,"setTXegAmpL1B",0));
  for(i=0;i<nWaves;i++)tempFloat[i]=maxAmp;

  return(tempFloat);
}/*setTXegAmpL1B*/


/*####################################################*/
/*set th_left_used for L1B files*/

uint16_t *setThUsedL1B(int nWaves)
{
  int i=0;
  uint16_t *tempUint16=NULL;

  if(!(tempUint16=(uint16_t *)calloc(nWaves,sizeof(uint16_t)))){
    errorf("error in tempUint16 allocation.\n");
    return NULL;
  }

  for(i=0;i<nWaves;i++)tempUint16[i]=260;

  return(tempUint16);
}/*setThUsedL1B*/


/*####################################################*/
/*an array of selection_stretchers_x*/

uint16_t *setSelectStretchL1B(int nWaves)
{
  int i=0;
  uint16_t *tempUint16=NULL;

  if(!(tempUint16=(uint16_t *)calloc(nWaves,sizeof(uint16_t)))){
    errorf("error in tempUint16 allocation.\n");
    return NULL;
  }

  for(i=0;i<nWaves;i++)tempUint16[i]=50;

  return(tempUint16);
}/*setSelectStretchL1B*/


/*####################################################*/
/*set RX energy*/

double *setRXenergyL1B(gediHDF *hdfData)
{
  int i=0,j=0;
  double *tempDouble=NULL;

  tempDouble=dalloc(hdfData->nWaves,"setRXenergyL1B",0);

  for(i=0;i<hdfData->nWaves;i++){
    tempDouble[i]=0.0;
    for(j=0;j<hdfData->nBins[0];j++)tempDouble[i]+=hdfData->wave[0][j];
  }

  return(tempDouble);
}/*setRXenergyL1B*/


/*####################################################*/
/*set all_samples_sum*/

uint32_t *setAllSamplesSumL1B(gediHDF *hdfData)
{
  int i=0,j=0;
  uint32_t *tempUint32=NULL;

  if(!(tempUint32=(uint32_t *)calloc(hdfData->nWaves,sizeof(uint32_t)))){
    fprintf(stderr,"error in tempUint32 allocation.\n");
    return(NULL);
  }

  for(i=0;i<hdfData->nWaves;i++){
    tempUint32[i]=0;
    for(j=0;j<hdfData->nBins[0];j++)tempUint32[i]+=hdfData->wave[0][j];
  }

  return(tempUint32);
}/*setAllSamplesSumL1B*/


/*####################################################*/
/*set neutat_delay_derivative*/

double *setDelayDerivL1B(int nWaves)
{
  int i=0;
  double *tempDouble=NULL;

  ASSIGN_CHECKNULL_RETNULL(tempDouble,dalloc(nWaves,"setDelayDerivL1B",0));
  for(i=0;i<nWaves;i++)tempDouble[i]=-8.82508928*pow(10.0,-8.0);

  return(tempDouble);
}/*setDelayDerivL1B*/


/*####################################################*/
/*set a 2D array of surface types for L1B format*/

int8_t *setSurfaceTypeL1B(int nWaves,int nLayers)
{
  int i=0,j=0;
  int8_t *tempInt8=NULL;

  /*allocate space*/
  if(!(tempInt8=(int8_t *)calloc(nWaves*nLayers,sizeof(int8_t)))){
    errorf("error in txCount allocation.\n");
    return NULL;
  } 

  for(i=0;i<nLayers;i++){
    if(i>0){
      for(j=0;j<nWaves;j++)tempInt8[i*nWaves+j]=1;
    }else{
      for(j=0;j<nWaves;j++)tempInt8[i*nWaves+j]=0;
    }
  }

  return(tempInt8);
}/*setSurfaceTypeL1B*/


/*####################################################*/
/*make an array of pi/2*/

double *setHalfPiL1B(int numb)
{
  int i=0;
  double *tempDouble=NULL,pi2=0;

  ASSIGN_CHECKNULL_RETNULL(tempDouble,dalloc(numb,"setHalfPiL1B",0));
  pi2=M_PI/2.0;
  for(i=0;i<numb;i++)tempDouble[i]=pi2;

  return(tempDouble);
}/*setHalfPiL1B*/


/*####################################################*/
/*reproject coordinates for L1B format*/
double *setL1Bcoords(int aEPSG,gediHDF *hdfData)
{
  #ifndef WITHOUT_GDAL
  int i=0;
  double *tempDouble=NULL,*z=NULL;
  double *reprojectWaveBounds(double *inBounds,int inEPSG,int outEPSG);
  OGRCoordinateTransformationH hTransform;
  OGRSpatialReferenceH hSourceSRS,hTargetSRS;

  /*allocate space*/
  ASSIGN_CHECKNULL_RETNULL(tempDouble,dalloc(2*hdfData->nWaves,"setL1Bcoords",0));  /*x*nWaves,y*nWaves*/

  /*copy relevant data*/
  for(i=0;i<hdfData->nWaves;i++){
    tempDouble[i]=hdfData->lon[i];
    tempDouble[i+hdfData->nWaves]=hdfData->lat[i];
  }

  /*does it need reprojecting?*/
  if(aEPSG!=4326){
    /*make dummy z array*/
    ASSIGN_CHECKNULL_RETNULL(z,dalloc(hdfData->nWaves,"Z setL1Bcoords",0));

    hSourceSRS=OSRNewSpatialReference(NULL);
    hTargetSRS=OSRNewSpatialReference(NULL);
    OSRImportFromEPSG(hTargetSRS,4326);
    OSRImportFromEPSG(hSourceSRS,aEPSG);
    hTransform=OCTNewCoordinateTransformation(hSourceSRS,hTargetSRS);
    OCTTransform(hTransform,hdfData->nWaves,&(tempDouble[0]),&(tempDouble[hdfData->nWaves]),z);
    OCTDestroyCoordinateTransformation(hTransform);
    OSRDestroySpatialReference(hSourceSRS);
    OSRDestroySpatialReference(hTargetSRS);
    TIDY(z);
  }

  return(tempDouble);
  #else
    return NULL;
  #endif
}/*setL1Bcoords*/

/*####################################################*/
/*set cover for L1B format*/

double *setCovL1B(gediHDF *hdfData,int useInt)
{
  int i=0,j=0;
  float tot=0,totG=0;
  double *tempDouble=NULL;

  ASSIGN_CHECKNULL_RETNULL(tempDouble,dalloc(hdfData->nWaves,"setCovL1B",0));

  if(hdfData->ground){
    /*loop over waves*/
    for(i=0;i<hdfData->nWaves;i++){
      tot=totG=0.0;
      for(j=0;j<hdfData->nBins[0];j++){
        tot+=hdfData->wave[useInt][j];
        totG+=hdfData->ground[useInt][j];
      }
      if(tot>0.0)tempDouble[i]=(double)((tot-totG)/tot)*100.0;  /*it is in pecent rather than a fraction*/
      else       tempDouble[i]=0.0;
    }
  }else for(i=0;i<hdfData->nWaves;i++)tempDouble[i]=50.0;       /*where no data, set as 50%*/

  return(tempDouble);
}/*setCovL1B*/


/*####################################################*/
/*set elevation bin for L1B format*/

double *setElevBinL1B(int nWaves,float *z0)
{
  int i=0;
  double *tempDouble=NULL;

  ASSIGN_CHECKNULL_RETNULL(tempDouble,dalloc(nWaves,"setElevBinL1B",0));
  for(i=0;i<nWaves;i++)tempDouble[i]=(double)z0[i];

  return(tempDouble);
}/*setElevBinL1B*/


/*####################################################*/
/*make a DEM for the l1b format*/

double *setDEMl1b(gediHDF *hdfData,int useInt)
{
  int i=0,j=0;
  double *tempDouble=NULL;
  double z=0,contN=0,res=0;

  /*allocate space*/
  ASSIGN_CHECKNULL_RETNULL(tempDouble,dalloc(hdfData->nWaves,"setDEMl1b",0));

  /*do we have ground data?*/
  if(hdfData->ground){
    /*loop over waves*/
    for(i=0;i<hdfData->nWaves;i++){
      tempDouble[i]=contN=0.0;
      res=(double)(hdfData->z0[i]-hdfData->zN[i])/(double)hdfData->nBins[0];
      /*calculate CofG*/
      for(j=0;j<hdfData->nBins[0];j++){
        z=(double)hdfData->z0[i]-(double)j*res;
        tempDouble[i]+=(double)hdfData->ground[useInt][j]*z;
        contN+=(double)hdfData->ground[useInt][j];
      }
      if(contN>0.0)tempDouble[i]/=contN;
    }
  }else for(i=0;i<hdfData->nWaves;i++)tempDouble[i]=0.0;

  return(tempDouble);
}/*setDEMl1b*/


/*####################################################*/
/*set delta time*/

double *setDeltaTime(int nWaves)
{
  int i=0;
  double *tempDouble=NULL;

  ASSIGN_CHECKNULL_RETNULL(tempDouble,dalloc(nWaves,"setDeltaTime",0));

  for(i=0;i<nWaves;i++)tempDouble[i]=40179291.0+(double)i/484.0;  /*a dummy time for now, but with the right spacing*/

  return(tempDouble);
}/*setDeltaTime*/


/*####################################################*/
/*set bounds offset for bins*/

double *setBounceOffset(int nWaves,int *nBins,float res)
{
  int i=0;
  double *tempDouble=NULL;
  double c=0;

  ASSIGN_CHECKNULL_RETNULL(tempDouble,dalloc(nWaves,"setBounceOffset",0));

  /*is this dfor bin 0 or bin 1*/
  if(nBins){
    c=299792458.0;
    for(i=0;i<nWaves;i++)tempDouble[i]=0.0013691+(double)nBins[0]*(double)res/c;
  }else{
    for(i=0;i<nWaves;i++)tempDouble[i]=0.0013691;
  }

  return(tempDouble);
}/*setBounceOffset*/


/*####################################################*/
/*set instrument altitude*/

double *setAltitude(int nWaves)
{
  int i=0;
  double *tempDouble=NULL;

  ASSIGN_CHECKNULL_RETNULL(tempDouble,dalloc(nWaves,"setAltitude",0));
  for(i=0;i<nWaves;i++)tempDouble[i]=410560.0;

  return(tempDouble);
}/*setAltitude*/


/*####################################################*/
/*Rearrange pulse information in to TX structure*/

int rearrangePulsetoTX(gediIOstruct *gediIO,gediHDF *hdfData,TXstruct *tx)
{
  int i=0,j=0;
  int *contN=NULL;
  int buff=0;
  int nPbins=0;
  float pRes=0,res=0;         /*pulse and wave resolution*/
  float *pulse=NULL,*pX=NULL; /*pointers to pulse amplitude and x*/
  float *setPx(float,int);
  float *txwave=NULL;         /*resamped pulse*/
  gediRatStruct gediRat;


  /*set pointers to piulse, depending on where it is defined*/
  if(gediIO->pulse){             /*pulse is defined in the sim settings*/
    pRes=gediIO->pRes;
    nPbins=gediIO->pulse->nBins;
    pulse=gediIO->pulse->y;
    pX=gediIO->pulse->x;
    res=gediIO->res;
  }else if(hdfData->pulse){      /*pulse is defined in a HDF file*/
    pRes=hdfData->pRes;
    nPbins=hdfData->nPbins;
    pulse=hdfData->pulse;
    pX=setPx(pRes,nPbins);
    res=(hdfData->z0[0]-hdfData->zN[0])/(float)(hdfData->nBins[0]-1);
  }else if(hdfData->pSigma>0.0){ /*pulse is Gaussian and needs defining*/
    /*set pulse*/
    gediIO->pSigma=hdfData->pSigma;
    gediIO->readPulse=0;
    gediRat.iThresh=0.0006;
    gediIO->res=(hdfData->z0[0]-hdfData->zN[0])/(float)(hdfData->nBins[0]-1);
    gediIO->pRes=gediIO->res;
    ISINTRETINT(setGediPulse(gediIO,&gediRat));
    /*set pointers*/
    pRes=gediIO->pRes;
    nPbins=gediIO->pulse->nBins;
    pulse=gediIO->pulse->y;
    pX=gediIO->pulse->x;
    res=gediIO->res;
  }else{
    errorf("No pulse defined. Cannot write L1B format without pulse.\n");
    return -1;
  }


  /*resample the resolution*/
  buff=50;   /*pad before and after the TX wave in case the L2A code needs some workspace*/
  tx->nBins=(uint16_t)((float)nPbins*pRes/res)+(uint16_t)(2*buff);

  /*allocate space*/
  if(!(tx->txCount=(uint16_t *)calloc(hdfData->nWaves,sizeof(uint16_t)))){
    errorf("error in txCount allocation.\n");
    return -1;
  }
  if(!(tx->txStart=(uint64_t *)calloc(hdfData->nWaves,sizeof(uint64_t)))){
    errorf("error in txStart allocation.\n");
    return -1;
  }


  /*make a resampled pulse*/
  ASSIGN_CHECKNULL_RETINT(tx->txwave,falloc((int)tx->nBins*hdfData->nWaves,"txwave",0));
  ASSIGN_CHECKNULL_RETINT(txwave,falloc((int)tx->nBins,"temp txwave",0));
  ASSIGN_CHECKNULL_RETINT(contN,ialloc((int)tx->nBins,"txwave counter",0));

  /*zero counters*/
  for(i=0;i<(int)tx->nBins;i++){
    txwave[i]=0.0;
    contN[i]=0;
  }


  /*count up*/

    for(i=0;i<nPbins;i++){
      j=(int)(pX[i]/res)+buff;
      if((j>=0)&&(j<(int)tx->nBins)){
        txwave[j]+=pulse[i];
        contN[j]++;
      }
    }

    /*normalise*/
    tx->maxAmp=0.0;
    for(i=0;i<(int)tx->nBins;i++){
      if(contN[i]>0){
        txwave[i]/=(float)contN[i];
        if(txwave[i]>tx->maxAmp)tx->maxAmp=txwave[i];
      }
    }

    /*populate arrays*/
    for(i=0;i<hdfData->nWaves;i++){
      tx->txCount[i]=tx->nBins;
      tx->txStart[i]=(uint64_t)i*(uint64_t)tx->nBins;
      memcpy(&(tx->txwave[i*(int)tx->nBins]),&(txwave[0]),sizeof(float)*tx->nBins);
    }
    TIDY(txwave);





  return 0;
}/*rearrangePulsetoTX*/


/*####################################################*/
/*set the x array for a pulse*/

float *setPx(float pRes,int nPbins)
{
  int i=0;
  float *pX=NULL;

  ASSIGN_CHECKNULL_RETNULL(pX,falloc(nPbins,"pulse x",0));
  for(i=0;i<nPbins;i++)pX[i]=(float)i*pRes;

  return(pX);
}/*setPx*/


/*####################################################*/
/*make an array of int32 of 1*/

int32_t *padInt32ones(int numb)
{
  int i=0;
  int32_t *jimlad=NULL;

  if(!(jimlad=(int32_t *)calloc(numb,sizeof(int32_t)))){
    errorf("error in padInt32ones allocation.\n");
    return NULL;
  }

  for(i=0;i<numb;i++)jimlad[i]=1;

  return(jimlad);
}/*padInt32ones*/


/*####################################################*/
/*set a shot number*/

uint64_t *setShotNumber(int nWaves)
{
  int i=0;
  uint64_t *tempUint64=NULL;

  if(!(tempUint64=(uint64_t *)calloc(nWaves,sizeof(uint64_t)))){
    errorf("error in tempUint64 allocation.\n");
    return NULL;
  }

  for(i=0;i<nWaves;i++)tempUint64[i]=(uint64_t)i;

  return(tempUint64);
}/*setShotNumber*/


/*####################################################*/
/*make array to map RX wave starts*/

uint64_t *setRXstarts(int nWaves,int *nBins)
{
  int i=0;
  uint64_t *tempUint64=NULL;

  if(!(tempUint64=(uint64_t *)calloc(nWaves,sizeof(uint64_t)))){
    errorf("error in tempUint64 allocation.\n");
    return NULL;
  }

  for(i=0;i<nWaves;i++)tempUint64[i]=(uint64_t)i*(uint64_t)nBins[0];

  return(tempUint64);
}/*setRXstarts*/


/*####################################################*/
/*get Rx Sample Count in to L1B format*/

uint16_t *setRxSampleCount(int *nBins,int nWaves)
{
  int i=0;
  uint16_t *tempUint16=NULL;

  if(!(tempUint16=(uint16_t *)calloc(nWaves,sizeof(uint16_t)))){
    errorf("error in tempUint16 allocation.\n");
    return NULL;
  }

  for(i=0;i<nWaves;i++)tempUint16[i]=(uint16_t)nBins[0];

  return(tempUint16);
}/*setRxSampleCount*/


/*####################################################*/
/*make an array of 0s in uint8 for HDF5*/

uint8_t *padUint8zeros(int numb)
{
  int i=0;
  uint8_t *tempUint8=NULL;

  if(!(tempUint8=(uint8_t *)calloc(numb,sizeof(uint8_t)))){
    errorf("error in tempUint8 allocation.\n");
    return NULL;
  } 

  for(i=0;i<numb;i++)tempUint8[i]=0;

  return(tempUint8);
}/*padUint8zeros*/


/*####################################################*/
/*make an array of 1s in uint8 for HDF5*/

uint8_t *padUint8ones(int numb)
{
  int i=0;
  uint8_t *tempUint8=NULL;

  if(!(tempUint8=(uint8_t *)calloc(numb,sizeof(uint8_t)))){
    errorf("error in tempUint8 allocation.\n");
    return NULL;
  }

  for(i=0;i<numb;i++)tempUint8[i]=1;

  return(tempUint8);
}/*padUint8ones*/


/*####################################################*/
/*make an array of 0s in uint16 for HDF5*/

uint16_t *padUint16zeros(int numb)
{
  int i=0;
  uint16_t *tempUint16=NULL;

  if(!(tempUint16=(uint16_t *)calloc(numb,sizeof(uint16_t)))){
    errorf("error in tempUint16 allocation.\n");
    return NULL;
  }

  for(i=0;i<numb;i++)tempUint16[i]=0;


  return(tempUint16);
}/*padUint16zeros*/


/*####################################################*/
/*make an array of 0s in uint32 for HDF5*/

uint32_t *padUint32zeros(int numb)
{
  int i=0;
  uint32_t *tempUint32=NULL;

  if(!(tempUint32=(uint32_t *)calloc(numb,sizeof(uint32_t)))){
    errorf("error in tempUint32 allocation.\n");
    return NULL;
  } 

  for(i=0;i<numb;i++)tempUint32[i]=0;

  return(tempUint32);
}/*padUint32zeros*/


/*####################################################*/
/*write data to HDF5*/

int writeGEDIhdf(gediHDF *hdfData,char *namen,gediIOstruct *gediIO)
{
  hid_t file;         /* Handles */

  /*open new file*/
  file=H5Fcreate(namen,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

 /*write header*/
  ISINTRETINT(write1dIntHDF5(file,"NWAVES",&hdfData->nWaves,1));
  ISINTRETINT(write1dIntHDF5(file,"NBINS",&hdfData->nBins[0],1));
  ISINTRETINT(write1dIntHDF5(file,"NTYPEWAVES",&hdfData->nTypeWaves,1));
  ISINTRETINT(write1dIntHDF5(file,"IDLENGTH",&hdfData->idLength,1));
  ISINTRETINT(write1dFloatHDF5(file,"FSIGMA",&hdfData->fSigma,1));
  ISINTRETINT(write1dFloatHDF5(file,"PSIGMA",&hdfData->pSigma,1));
  ISINTRETINT(write1dIntHDF5(file,"NPBINS",&hdfData->nPbins,1));
  /*write datasets*/
  ISINTRETINT(write1dDoubleHDF5(file,"LON0",hdfData->lon,hdfData->nWaves));
  ISINTRETINT(write1dDoubleHDF5(file,"LAT0",hdfData->lat,hdfData->nWaves));
  if(hdfData->ground){
    ISINTRETINT(writeComp1dFloatHDF5(file,"SLOPE",hdfData->slope,hdfData->nWaves));
    ISINTRETINT(writeComp1dFloatHDF5(file,"ZG",hdfData->gElev,hdfData->nWaves));
    if(hdfData->demElev)writeComp1dFloatHDF5(file,"ZGDEM",hdfData->demElev,hdfData->nWaves);
  }
  ISINTRETINT(writeComp1dFloatHDF5(file,"BEAMDENSE",hdfData->beamDense,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(file,"POINTDENSE",hdfData->pointDense,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(file,"INCIDENTANGLE",hdfData->zen,hdfData->nWaves));
  if(gediIO->useInt){
    ISINTRETINT(writeComp2dFloatHDF5(file,"RXWAVEINT",hdfData->wave[0],hdfData->nWaves,hdfData->nBins[0]));
    if(hdfData->ground) {ISINTRETINT(writeComp2dFloatHDF5(file,"GRWAVEINT",hdfData->ground[0],hdfData->nWaves,hdfData->nBins[0]));}
  }
  if(gediIO->useCount){
    ISINTRETINT(writeComp2dFloatHDF5(file,"RXWAVECOUNT",hdfData->wave[(int)gediIO->useInt],hdfData->nWaves,hdfData->nBins[0]));
    if(hdfData->ground) {ISINTRETINT(writeComp2dFloatHDF5(file,"GRWAVECOUNT",hdfData->ground[(int)gediIO->useInt],hdfData->nWaves,hdfData->nBins[0]));}
  }
  if(gediIO->useFrac){
    ISINTRETINT(writeComp2dFloatHDF5(file,"RXWAVEFRAC",hdfData->wave[(int)(gediIO->useCount+gediIO->useFrac)],hdfData->nWaves,hdfData->nBins[0]));
    if(hdfData->ground) {ISINTRETINT(writeComp2dFloatHDF5(file,"GRWAVEFRAC",hdfData->ground[(int)(gediIO->useCount+gediIO->useFrac)],hdfData->nWaves,hdfData->nBins[0]));}
  }
  ISINTRETINT(writeComp1dFloatHDF5(file,"Z0",hdfData->z0,hdfData->nWaves));
  ISINTRETINT(writeComp1dFloatHDF5(file,"ZN",hdfData->zN,hdfData->nWaves));
  ISINTRETINT(write2dCharHDF5(file,"WAVEID",hdfData->waveID,hdfData->nWaves,hdfData->idLength));

  if(hdfData->nPbins>0){
    ISINTRETINT(write1dFloatHDF5(file,"PRES",&hdfData->pRes,1));
    ISINTRETINT(writeComp1dFloatHDF5(file,"PULSE",hdfData->pulse,hdfData->nPbins));
  }

  /*close file*/
  if(H5Fclose(file)){
    errorf("Issue closing file\n");
    return(-1);
  }
  msgf("Waveforms written to %s\n",namen);
  return(0);
}/*writeGEDIhdf*/


/*####################################################*/
/*read HDF5 GEDI data into structure*/

gediHDF *readGediHDF(char *namen,gediIOstruct *gediIO)
{
  hid_t file;         /* Handles */
  gediHDF *hdfData=NULL;
  int readSimGediHDF(hid_t,gediIOstruct *,char *,gediHDF *);
  int readRealGediHDF(hid_t,gediIOstruct *,char *,gediHDF *);

  /*allocate space for all*/
  if(!(hdfData=(gediHDF *)calloc(1,sizeof(gediHDF)))){
    errorf("error control allocation.\n");
    return(NULL);
  }

  /*open file*/
  msgf("Reading %s\n",namen);
  file=H5Fopen(namen,H5F_ACC_RDONLY,H5P_DEFAULT);

  /*is the file a simulation or real?*/
  int result=0;
  if(H5Lexists(file,"BEAMDENSE",H5P_DEFAULT)==0)result=readRealGediHDF(file,gediIO,namen,hdfData);
  else                                          result=readSimGediHDF(file,gediIO,namen,hdfData);
  if(result==-1)
    return(NULL);

  /*close file*/
  if(H5Fclose(file)){
    errorf("Issue closing file\n");
    return(NULL);
  }
  return(hdfData);
}/*readGediHDF*/


/*####################################################*/
/*read simulated HDF GEDI data*/

int readSimGediHDF(hid_t file,gediIOstruct *gediIO,char *namen,gediHDF *hdfData)
{
  int i=0;
  int nWaves=0,nBins=0;
  int *tempI=NULL,ind=0;
  float *tempF=NULL;

  /*sims do not have variable lengths*/
  ASSIGN_CHECKNULL_RETINT(hdfData->nBins,ialloc(1,"nBins",0));
  hdfData->varBins=0;

  /*read the header*/
  ASSIGN_CHECKNULL_RETINT(tempF,read1dFloatHDF5(file,"FSIGMA",&nWaves));
  hdfData->fSigma=*tempF;
  TIDY(tempF);
  ASSIGN_CHECKNULL_RETINT(tempF,read1dFloatHDF5(file,"PSIGMA",&nWaves));
  hdfData->pSigma=*tempF;
  TIDY(tempF);
  ASSIGN_CHECKNULL_RETINT(tempI,read1dIntHDF5(file,"NWAVES",&nWaves));
  hdfData->nWaves=*tempI;
  TIDY(tempI);
  ASSIGN_CHECKNULL_RETINT(tempI,read1dIntHDF5(file,"NBINS",&nWaves));
  hdfData->nBins[0]=*tempI;
  TIDY(tempI);
  ASSIGN_CHECKNULL_RETINT(tempI,read1dIntHDF5(file,"NPBINS",&nWaves));
  hdfData->nPbins=*tempI;
  TIDY(tempI);
  ASSIGN_CHECKNULL_RETINT(tempI,read1dIntHDF5(file,"NTYPEWAVES",&nWaves));
  hdfData->nTypeWaves=*tempI;
  TIDY(tempI);
  ASSIGN_CHECKNULL_RETINT(tempI,read1dIntHDF5(file,"IDLENGTH",&nWaves));
  hdfData->idLength=*tempI;
  TIDY(tempI);

  /*set which waves are here if not predefined*/
  if(!(gediIO->useInt+gediIO->useCount+gediIO->useFrac)){
    if(hdfData->nTypeWaves==3){
      gediIO->useInt=gediIO->useCount=gediIO->useFrac=1;
    }else{
      gediIO->useInt=gediIO->useFrac=0;
      gediIO->useCount=1;
    }
  }

  /*read ancillary data*/
  ASSIGN_CHECKNULL_RETINT(hdfData->z0,read1dFloatHDF5(file,"Z0",&nWaves));
  ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"Z0"));
  ASSIGN_CHECKNULL_RETINT(hdfData->zN,read1dFloatHDF5(file,"ZN",&nWaves));
  ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"ZN"));
  if(gediIO->ground){
    ASSIGN_CHECKNULL_RETINT(hdfData->slope,read1dFloatHDF5(file,"SLOPE",&nWaves));
    ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"SLOPE"));
    ASSIGN_CHECKNULL_RETINT(hdfData->gElev,read1dFloatHDF5(file,"ZG",&nWaves));
    ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"ZG"));
  }
  ASSIGN_CHECKNULL_RETINT(hdfData->beamDense,read1dFloatHDF5(file,"BEAMDENSE",&nWaves));
  ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"BEAMDENSE"));
  ASSIGN_CHECKNULL_RETINT(hdfData->pointDense,read1dFloatHDF5(file,"POINTDENSE",&nWaves));
  ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"POINTDENSE"));
  ASSIGN_CHECKNULL_RETINT(hdfData->zen,read1dFloatHDF5(file,"INCIDENTANGLE",&nWaves));
  ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"INCIDENTANGLE"));
  ASSIGN_CHECKNULL_RETINT(hdfData->lon,read1dDoubleHDF5(file,"LON0",&nWaves));
  ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"LON0"));
  ASSIGN_CHECKNULL_RETINT(hdfData->lat,read1dDoubleHDF5(file,"LAT0",&nWaves));
  ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"LAT0"));
  ASSIGN_CHECKNULL_RETINT(hdfData->waveID,read15dCharHDF5(file,"WAVEID",&nWaves,&nBins));
  ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"WAVEID"));
  if(hdfData->nPbins>0){
    ASSIGN_CHECKNULL_RETINT(hdfData->pulse,read1dFloatHDF5(file,"PULSE",&nBins));
    ISINTRETINT(checkNwavesDF(nBins,hdfData->nPbins,"PULSE"));
    ASSIGN_CHECKNULL_RETINT(tempF,read1dFloatHDF5(file,"PRES",&nWaves));
    hdfData->pRes=*tempF;
    TIDY(tempF);
  }else hdfData->pulse=NULL;

  //float *demElev;  /*ground elevation, DEM*/


  /*determine how many waveforms we want to read*/
  if((int)(gediIO->useInt+gediIO->useCount+gediIO->useFrac)>hdfData->nTypeWaves){
    errorf("Not enough waveform types for that option set. %d %d from %d %d %d\n",hdfData->nTypeWaves,\
               gediIO->useInt+gediIO->useCount+gediIO->useFrac,gediIO->useCount,gediIO->useInt,gediIO->useFrac);
    return(-1);
  }else hdfData->nTypeWaves=(int)(gediIO->useInt+gediIO->useCount+gediIO->useFrac);

  /*allocate waveform space*/
  ASSIGN_CHECKNULL_RETINT(hdfData->wave,fFalloc(hdfData->nTypeWaves,"hdf waveforms",0));
  if(gediIO->ground) {ASSIGN_CHECKNULL_RETINT(hdfData->ground,fFalloc(hdfData->nTypeWaves,"hdf waveforms",0));}
  if(!(hdfData->sInd=(uint64_t *)calloc(hdfData->nWaves,sizeof(uint64_t)))){
    errorf("error in sInd buffer allocation.\n");
    return(-1);
  }
  for(i=0;i<hdfData->nWaves;i++)hdfData->sInd[i]=(uint64_t)i*(uint64_t)hdfData->nBins[0];

  /*read data*/
  if(gediIO->useInt){
    ind=0;
    ASSIGN_CHECKNULL_RETINT(hdfData->wave[ind],read15dFloatHDF5(file,"RXWAVEINT",&nWaves,&nBins));
    ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"RXWAVEINT"));
    ISINTRETINT(checkNwavesDF(nBins,hdfData->nBins[0],"RXWAVEINT"));
    if(gediIO->ground){
      ASSIGN_CHECKNULL_RETINT(hdfData->ground[ind],read15dFloatHDF5(file,"GRWAVEINT",&nWaves,&nBins));
      ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"GRWAVEINT"));
      ISINTRETINT(checkNwavesDF(nBins,hdfData->nBins[0],"GRWAVEINT"));
    }
  }
  if(gediIO->useCount){
    ind=(int)gediIO->useInt;
    ASSIGN_CHECKNULL_RETINT(hdfData->wave[ind],read15dFloatHDF5(file,"RXWAVECOUNT",&nWaves,&nBins));
    ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"RXWAVECOUNT"));
    ISINTRETINT(checkNwavesDF(nBins,hdfData->nBins[0],"RXWAVECOUNT"));
    if(gediIO->ground){
      ASSIGN_CHECKNULL_RETINT(hdfData->ground[ind],read15dFloatHDF5(file,"GRWAVECOUNT",&nWaves,&nBins));
      ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"GRWAVECOUNT"));
      ISINTRETINT(checkNwavesDF(nBins,hdfData->nBins[0],"GRWAVECOUNT"));
    }
  }
  if(gediIO->useFrac){
    ind=(int)(gediIO->useInt+gediIO->useCount);
    ASSIGN_CHECKNULL_RETINT(hdfData->wave[ind],read15dFloatHDF5(file,"RXWAVEFRAC",&nWaves,&nBins));
    ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"RXWAVEFRAC"));
    ISINTRETINT(checkNwavesDF(nBins,hdfData->nBins[0],"RXWAVEFRAC"));
    if(gediIO->ground){
      ASSIGN_CHECKNULL_RETINT(hdfData->ground[ind],read15dFloatHDF5(file,"GRWAVEFRAC",&nWaves,&nBins));
      ISINTRETINT(checkNwavesDF(nWaves,hdfData->nWaves,"GRWAVEFRAC"));
      ISINTRETINT(checkNwavesDF(nBins,hdfData->nBins[0],"GRWAVEFRAC"));
    }
  }

  /*some parameters are not in the simulations*/
  hdfData->solarElev=NULL;

  return(0);
}/*readSimGediHDF*/


/*####################################################*/
/*read real HDF GEDI data*/

int readRealGediHDF(hid_t file,gediIOstruct *gediIO,char *namen,gediHDF *hdfData)
{
  int i=0,j=0,nBeams=0;
  int numb=0,nSamps=0;
  int *useInd=NULL,nUse=0;
  int *usableGEDIfootprints(double *,double *,int,int *,gediIOstruct *);
  uint16_t *nBins=NULL;
  uint64_t *sInds=NULL;
  uint64_t *shotN=NULL;
  hid_t group=0,group2=0;
  double *temp1=NULL,*temp2=NULL;
  double *tempLon=NULL,*tempLat=NULL;
  double *meanCoord(double *,double *,int);
  char l1b=0;     /*l1a or l1b flag*/
  char **beamList=NULL;
  char **setGEDIbeamList(int *,char *);
  int updateGEDInWaves(int,gediHDF *);
  void setGEDIzenith(gediHDF *,int,uint16_t *);
  void setGEDIwaveID(gediHDF *,int,uint64_t *,int *,char *);
  int readGEDIwaveform(hid_t,int *,uint64_t *,int,gediHDF *,int *,char *);

  ASSIGN_CHECKNULL_RETINT(beamList,setGEDIbeamList(&nBeams,gediIO->useBeam));
  /*set the list of beams*/
  hdfData->nWaves=0;
  hdfData->varBins=1;  /*variable bin length*/
  gediIO->ground=0;    /*no truth with real data*/
  gediIO->nTypeWaves=hdfData->nTypeWaves=1;
  gediIO->useCount=1;
  gediIO->useFrac=gediIO->useInt=0;


  /*loop over beams and read all*/
  for(i=0;i<nBeams;i++){

    /*does this beam exist in this file?*/
    if(H5Lexists(file,beamList[i],H5P_DEFAULT)==0)continue;

    /*open beam group*/
    group=H5Gopen2(file,beamList[i],H5P_DEFAULT);

    /*check whether this beam has data*/
    if(H5Lexists(group,"geolocation",H5P_DEFAULT)==0)continue;

    /*geolocation*/
    group2=H5Gopen2(group,"geolocation",H5P_DEFAULT);
    ASSIGN_CHECKNULL_RETINT(temp1,read1dDoubleHDF5(group2,"longitude_bin0",&numb));
    ASSIGN_CHECKNULL_RETINT(temp2,read1dDoubleHDF5(group2,"longitude_lastbin",&numb));
    ASSIGN_CHECKNULL_RETINT(tempLon,meanCoord(temp1,temp2,numb));
    TIDY(temp1);
    TIDY(temp2);
    ASSIGN_CHECKNULL_RETINT(temp1,read1dDoubleHDF5(group2,"latitude_bin0",&numb));
    ASSIGN_CHECKNULL_RETINT(temp2,read1dDoubleHDF5(group2,"latitude_lastbin",&numb));
    ASSIGN_CHECKNULL_RETINT(tempLat,meanCoord(temp1,temp2,numb));
    TIDY(temp1);
    TIDY(temp2);


    /*which are within bounds?*/
    ASSIGN_CHECKNULL_RETINT(useInd,usableGEDIfootprints(tempLon,tempLat,numb,&nUse,gediIO));

    /*are there any usable in this track?*/
    if(nUse>0){
      /*update the number of waves and arrays for holding data*/
      ISINTRETINT(updateGEDInWaves(nUse,hdfData));

      for(j=0;j<nUse;j++){
        hdfData->lat[j+hdfData->nWaves]=tempLat[useInd[j]];
        hdfData->lon[j+hdfData->nWaves]=tempLon[useInd[j]];
      }
      TIDY(tempLon);
      TIDY(tempLat);

      ASSIGN_CHECKNULL_RETINT(temp1,read1dDoubleHDF5(group2,"elevation_bin0",&numb));
      for(j=0;j<nUse;j++)hdfData->z0[j+hdfData->nWaves]=(float)temp1[useInd[j]];
      TIDY(temp1);
      ASSIGN_CHECKNULL_RETINT(temp1,read1dDoubleHDF5(group2,"elevation_lastbin",&numb));
      for(j=0;j<nUse;j++)hdfData->zN[j+hdfData->nWaves]=(float)temp1[useInd[j]];
      TIDY(temp1);

      /*waveform*/
      ASSIGN_CHECKNULL_RETINT(nBins,read1dUint16HDF5(group,"rx_sample_count",&numb));
      for(j=0;j<nUse;j++)hdfData->nBins[j+hdfData->nWaves]=(int)nBins[useInd[j]];
      TIDY(nBins);
      ASSIGN_CHECKNULL_RETINT(sInds,read1dUint64HDF5(group,"rx_sample_start_index",&numb));
      ISINTRETINT(readGEDIwaveform(group,&nSamps,sInds,nUse,hdfData,useInd,&l1b));
      TIDY(sInds);

      /*calculate zenith angles from elevations*/
      setGEDIzenith(hdfData,numb,nBins);
      TIDY(nBins);

      /*if it is an L1B file read the solar elevation*/
      if(l1b){
        ASSIGN_CHECKNULL_RETINT(temp1,read1dDoubleHDF5(group2,"solar_elevation",&numb));
        for(j=0;j<nUse;j++)hdfData->solarElev[j+hdfData->nWaves]=(float)temp1[useInd[j]];
        TIDY(temp1);
      }
      H5Gclose(group2);

      /*set waveIDs*/
      ASSIGN_CHECKNULL_RETINT(shotN,read1dUint64HDF5(group,"shot_number",&numb));
      setGEDIwaveID(hdfData,nUse,shotN,useInd,beamList[i]);
      TIDY(shotN);

      hdfData->nWaves+=nUse;
    }else H5Gclose(group2);
    H5Gclose(group);
    TIDY(useInd);
  }/*beam loop*/

  TTIDY((void **)beamList,nBeams);

  if(hdfData->nWaves==0){
    errorf("No footprints contained\n");
    return(-1);
  }
  return(0);
}/*readRealGediHDF*/


/*####################################################*/
/*read a real GEDI waveform*/

int readGEDIwaveform(hid_t group,int *nSamps,uint64_t *sInds,int nUse,gediHDF *hdfData,int *useInd,char *l1b)
{
  uint16_t *tempI=NULL;
  float *tempF=NULL;
  hid_t dset,dtype;
  herr_t status;
  int unwrapRealGEDI(uint16_t *,float *,uint64_t *,int,int,gediHDF *,int *);

  /*read data dtype*/
  dset=H5Dopen2(group,"rxwaveform",H5P_DEFAULT);
  dtype=H5Dget_type(dset);
  status=H5Dclose(dset);
 if(status<0){
    errorf("Error closing HDF5 group\n");
    return -1;
  }

  /*checkdata type and read in to appropriate array*/
  if(H5Tequal(dtype,H5T_NATIVE_USHORT)||H5Tequal(dtype,H5T_NATIVE_UINT16)){   /*l1a file*/
    ASSIGN_CHECKNULL_RETINT(tempI,read1dUint16HDF5(group,"rxwaveform",nSamps));
    *l1b=0;
  }else if(H5Tequal(dtype,H5T_NATIVE_FLOAT)||H5Tequal(dtype,H5T_STD_I32BE)||H5Tequal(dtype,H5T_STD_I32LE)||\
            H5Tequal(dtype,H5T_IEEE_F32BE)||H5Tequal(dtype,H5T_IEEE_F32LE)||H5Tequal(dtype,H5T_INTEL_I32)||H5Tequal(dtype,H5T_INTEL_B32)){   /*l1b file*/
    ASSIGN_CHECKNULL_RETINT(tempF,read1dFloatHDF5(group,"rxwaveform",nSamps));
    *l1b=1;
  }else{
    errorf("rxwaveform data type not recognised\n");
    return(-1);
  }
  status=H5Tclose(dtype);
 if(status<0){
    errorf("Error closing HDF5 group\n");
    return -1;
  }

  /*unpack and pad all waves to have the same number of bins*/
  status=unwrapRealGEDI(tempI,tempF,sInds,*nSamps,nUse,hdfData,useInd);
  if(status==-1)
    return(-1);

  TIDY(tempI);
  TIDY(tempF);
  return(0);
}/*readGEDIwaveform*/


/*####################################################*/
/*set waveform ID*/

void setGEDIwaveID(gediHDF *hdfData,int nUse,uint64_t *shotN,int *useInd,char *beam)
{
  int i=0,ind=0;

  for(i=0;i<nUse;i++){
    ind=i+hdfData->nWaves;
    snprintf(&hdfData->waveID[ind*hdfData->idLength], hdfData->idLength,"gedi.%s.%" PRIu64,beam,shotN[useInd[i]]);
  }

  return;
}/*setGEDIwaveID*/


/*####################################################*/
/*determine which are within bounds*/

int *usableGEDIfootprints(double *tempLon,double *tempLat,int numb,int *nUse,gediIOstruct *gediIO)
{
  int i=0;
  int *useInd=NULL;
  int *markInt(int,int *,int);
  double *bounds=NULL;

  /*reset counter*/
  *nUse=0;

  /*reproject bounds if needed*/
  ASSIGN_CHECKNULL_RETNULL(bounds,reprojectWaveBounds(&(gediIO->bounds[0]),gediIO->bEPSG,gediIO->wEPSG));

  /*are bounds being used?*/
  if(fabs(bounds[0]-bounds[1])<TOL){
    bounds[0]=bounds[1]=-100000000.0;
    bounds[2]=bounds[3]=100000000.0;
  }

  /*if bounds are in degrees, do we need to unwrap?*/
  if(gediIO->wEPSG==4326){
    if(bounds[0]<0.0)bounds[0]+=360.0;
    if(bounds[2]<0.0)bounds[2]+=360.0;
  }

  /*loop over all footprints*/
  for(i=0;i<numb;i++){
    /*may need to unwrap longitudes here too*/
    if(gediIO->wEPSG==4326)if(tempLon[i]<0.0)tempLon[i]+=360.0;

    if((tempLon[i]>=bounds[0])&&(tempLon[i]<=bounds[2])&&\
       (tempLat[i]>=bounds[1])&&(tempLat[i]<=bounds[3])){
      ASSIGN_CHECKNULL_RETNULL(useInd,markInt(*nUse,useInd,i));
      (*nUse)++;
    }
  }


  TIDY(bounds);
  return(useInd);
}/*usableGEDIfootprints*/


/*####################################################*/
/*reproject waveform bounds if needed*/

double *reprojectWaveBounds(double *inBounds,int inEPSG,int outEPSG)
{
#ifndef WITHOUT_GDAL
  double *x=NULL,*y=NULL,*z=NULL;
  OGRCoordinateTransformationH hTransform;
  OGRSpatialReferenceH hSourceSRS,hTargetSRS;
  double *bounds=NULL;
  int verMaj=0;
  int findGDAlVerMaj();


  /*allocate space*/
  ASSIGN_CHECKNULL_RETNULL(bounds,dalloc(4,"wave bounds",0));

  /*do we need to transform?*/
  if(inEPSG!=outEPSG){
    ASSIGN_CHECKNULL_RETNULL(x,dalloc(2,"x trans",9));
    ASSIGN_CHECKNULL_RETNULL(y,dalloc(2,"y trans",9));
    ASSIGN_CHECKNULL_RETNULL(z,dalloc(2,"z trans",9));

    x[0]=inBounds[0];
    y[0]=inBounds[1];
    x[1]=inBounds[2];
    y[1]=inBounds[3];
    z[0]=z[1]=0.0;

    hSourceSRS=OSRNewSpatialReference(NULL);
    hTargetSRS=OSRNewSpatialReference(NULL);
    OSRImportFromEPSG(hTargetSRS,outEPSG);
    OSRImportFromEPSG(hSourceSRS,inEPSG);
    hTransform=OCTNewCoordinateTransformation(hSourceSRS,hTargetSRS);
    OCTTransform(hTransform,2,x,y,z);
    OCTDestroyCoordinateTransformation(hTransform);
    OSRDestroySpatialReference(hSourceSRS);
    OSRDestroySpatialReference(hTargetSRS);

    /*GDAL 3.0 and later now returns lat lon rather than lon lat. Find majer version*/
    /*this will need updating once we hit version 10*/
    verMaj=findGDAlVerMaj();

    if(verMaj>=3){  /*if GDAL >=v3, need to swap lat and lon*/
      bounds[0]=y[0];
      bounds[1]=x[0];
      bounds[2]=y[1];
      bounds[3]=x[1];
    }else{
      bounds[0]=x[0];
      bounds[1]=y[0];
      bounds[2]=x[1];
      bounds[3]=y[1];
    }
  }else{  /*copy bounds*/
    bounds[0]=inBounds[0];
    bounds[1]=inBounds[1];
    bounds[2]=inBounds[2];
    bounds[3]=inBounds[3];
  }

  TIDY(x);
  TIDY(y);
  TIDY(z);
  return(bounds);
#else
  double* bounds = NULL;
  ASSIGN_CHECKNULL_RETNULL(bounds, dalloc(4, "wave bounds", 0));
  bounds[0] = inBounds[0];
  bounds[1] = inBounds[1];
  bounds[2] = inBounds[2];
  bounds[3] = inBounds[3];
  return(bounds);
#endif
}/*reprojectWaveBounds*/


/*####################################################*/
/*find GDAL version major*/

int findGDAlVerMaj()
{
  int verMaj=0;
  #ifndef WITHOUT_GDAL
  float val=0;
  char vers[20];   /*GDAL version number string*/

  strcpy(&(vers[0]),GDALVersionInfo("VERSION_NUM"));
  val=atof(vers);
  verMaj=(int)(val/pow(10,(int)(log(val)/log(10.0))));

  #endif
  return(verMaj);
}/*findGDAlVerMaj*/


/*####################################################*/
/*unwrap real GEDI data*/

int unwrapRealGEDI(uint16_t *tempI,float *tempF,uint64_t *sInds,int nSamps,int nUse,gediHDF *hdfData,int *useInd)
{
  int i=0,j=0,ind=0;
  uint64_t totBins=0;
  uint64_t tPlace=0;
  uint64_t offset=0;

  /*number of bins in this section and find offsets*/
  totBins=0;
  for(i=0;i<nUse;i++){
    ind=i+hdfData->nWaves;
    totBins+=(uint64_t)hdfData->nBins[ind];
    if(ind>0)hdfData->sInd[ind]=hdfData->sInd[ind-1]+(uint64_t)hdfData->nBins[ind-1];
    else     hdfData->sInd[ind]=0;
  }

  /*allocate space*/
  if(hdfData->nWaves==0){
    ASSIGN_CHECKNULL_RETINT(hdfData->wave,fFalloc(1,"wave",0));
    ASSIGN_CHECKNULL_RETINT(hdfData->wave[0],falloc((uint64_t)totBins,"wave",0));
    offset=0;
  }else{
    offset=hdfData->sInd[hdfData->nWaves];
    if(!(hdfData->wave[0]=(float *)realloc(hdfData->wave[0],(totBins+offset)*(uint64_t)sizeof(float)))){
      errorf("Error in reallocation, allocating %" PRIu64 "\n",(totBins+offset)*(uint64_t)sizeof(float *));
      return(-1);
    }
  }

  /*copy data*/
  if(tempI){   /*reading uint16 or float?*/
    for(i=0;i<nUse;i++){
      ind=i+hdfData->nWaves;
      for(j=0;j<hdfData->nBins[ind];j++){
        tPlace=sInds[useInd[i]]+(uint64_t)j;
        hdfData->wave[0][offset]=(float)tempI[tPlace];
        offset++;
      }
    }
  }else{
    for(i=0;i<nUse;i++){
      ind=i+hdfData->nWaves;
      for(j=0;j<hdfData->nBins[ind];j++){
        tPlace=sInds[useInd[i]]+(uint64_t)j;
        hdfData->wave[0][offset]=tempF[tPlace];
        offset++;
      }
    }
  }

  return(0);
}/*unwrapRealGEDI*/


/*####################################################*/
/*calculate zenith angle for real GEDI data*/

void setGEDIzenith(gediHDF *hdfData,int numb,uint16_t *nBins)
{
  /* To be finished
  int i=0,ind=0;

  for(i=0;i<numb;i++){
    ind=i-hdfData->nWaves;

  }*/

  return;
}/*setGEDIzenith*/


/*####################################################*/
/*update the number of GEDI waves*/

int updateGEDInWaves(int numb,gediHDF *hdfData)
{
  /*if already allocated, reallocate*/
  if(hdfData->nWaves>0){
    if(!(hdfData->z0=(float *)realloc(hdfData->z0,((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(float)))){
      errorf("Error in reallocation, allocating %" PRIu64 "\n",((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(float));
      return(-1);
    }
    if(!(hdfData->zN=(float *)realloc(hdfData->zN,((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(float)))){
      errorf("Error in reallocation, allocating %" PRIu64 "\n",((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(float));
      return(-1);
    }
    if(!(hdfData->lon=(double *)realloc(hdfData->lon,((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(double)))){
      errorf("Error in reallocation, allocating %" PRIu64 "\n",((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(double));
      return(-1);
    }
    if(!(hdfData->lat=(double *)realloc(hdfData->lat,((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(double)))){
      errorf("Error in reallocation, allocating %" PRIu64 "\n",((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(double));
      return(-1);
    }
    if(!(hdfData->zen=(float *)realloc(hdfData->zen,((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(float)))){
      errorf("Error in reallocation, allocating %" PRIu64 "\n",((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(float));
      return(-1);
    }
    if(!(hdfData->solarElev=(float *)realloc(hdfData->solarElev,((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(float)))){
      errorf("Error in reallocation, allocating %" PRIu64 "\n",((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(float));
      return(-1);
    }
    if(!(hdfData->nBins=(int *)realloc(hdfData->nBins,((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(int)))){
      errorf("Error in reallocation, allocating %" PRIu64 "\n",((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(float));
      return(-1);
    }
    if(!(hdfData->sInd=(uint64_t *)realloc(hdfData->sInd,((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(uint64_t)))){
      errorf("Error in reallocation, allocating %" PRIu64 "\n",((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(float));
      return(-1);
    }
    if(!(hdfData->waveID=(char *)realloc(hdfData->waveID,((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)hdfData->idLength*(uint64_t)sizeof(char)))){
      errorf("Error in reallocation, allocating %" PRIu64 "\n",((uint64_t)numb+(uint64_t)hdfData->nWaves)*(uint64_t)sizeof(float));
      return(-1);
    }
  }else{  /*allocate for the first time*/
    ASSIGN_CHECKNULL_RETINT(hdfData->z0,falloc((uint64_t)numb,"z0",0));
    ASSIGN_CHECKNULL_RETINT(hdfData->zN,falloc((uint64_t)numb,"zN",0));
    ASSIGN_CHECKNULL_RETINT(hdfData->lon,dalloc(numb,"lon",0));
    ASSIGN_CHECKNULL_RETINT(hdfData->lat,dalloc(numb,"lat",0));
    ASSIGN_CHECKNULL_RETINT(hdfData->waveID,challoc(numb,"waveID",0));
    ASSIGN_CHECKNULL_RETINT(hdfData->zen,falloc((uint64_t)numb,"zen",0));
    ASSIGN_CHECKNULL_RETINT(hdfData->solarElev,falloc((uint64_t)numb,"solarElev",0));
    ASSIGN_CHECKNULL_RETINT(hdfData->nBins,ialloc(numb,"nBins",0));
    if(!(hdfData->sInd=(uint64_t *)calloc(numb,sizeof(uint64_t)))){
      errorf("error in sInd allocation.\n");
      return(-1);
    }
    hdfData->idLength=50;
    ASSIGN_CHECKNULL_RETINT(hdfData->waveID,challoc(numb*hdfData->idLength,"waveID",0));
    hdfData->ground=NULL;
    hdfData->slope=NULL;
    hdfData->gElev=NULL;
    hdfData->demElev=NULL;
    hdfData->beamDense=NULL;
    hdfData->pointDense=NULL;
    hdfData->nTypeWaves=1;
  }
  return(0);
}/*updateGEDInWaves*/


/*####################################################*/
/*get mean coordinate*/

double *meanCoord(double *temp1,double *temp2,int numb)
{
  int i=0;
  double *coord=NULL;

  ASSIGN_CHECKNULL_RETNULL(coord,dalloc(numb,"mean coord",0));
  for(i=0;i<numb;i++)coord[i]=(temp1[i]+temp2[i])/2.0;

  return(coord);
}/*meanCoord*/


/*####################################################*/
/*set list of GEDI beams*/

char **setGEDIbeamList(int *nBeams,char *useBeam)
{
  int i=0,count=0;
  char **beamList=NULL;
  char tempBeam[8][9];

  /*temporary list of all beams*/
  strcpy(tempBeam[0],"BEAM0000");
  strcpy(tempBeam[1],"BEAM0001");
  strcpy(tempBeam[2],"BEAM0010");
  strcpy(tempBeam[3],"BEAM0011");
  strcpy(tempBeam[4],"BEAM0101");
  strcpy(tempBeam[5],"BEAM0110");
  strcpy(tempBeam[6],"BEAM1000");
  strcpy(tempBeam[7],"BEAM1011");

  /*how many beams shall we use?*/
  *nBeams=0;
  for(i=0;i<8;i++)if(useBeam[i])(*nBeams)++;

  /*allocate*/
  ASSIGN_CHECKNULL_RETNULL(beamList,chChalloc(*nBeams,"beam list",0));
  for(i=0;i<(*nBeams);i++) {ASSIGN_CHECKNULL_RETNULL(beamList[i],challoc(9,"beam list",i+1));}

  /*copy over*/
  count=0;
  for(i=0;i<8;i++){
    if(useBeam[i]){
      strcpy(beamList[count],tempBeam[i]);
      count++;
    }
  }

  return(beamList);
}/*setGEDIbeamList*/


/*####################################################*/
/*read an 8 bit strong and turn into an array*/

void setBeamsToUse(char *useBeam,char *instruction)
{ 
  int i=0,count=0;
  char temp[1];

  /*count down for efficiency*/
  count=7;
  for(i=strlen(instruction)-1;i>=0;i--){
    if((!strncasecmp(&(instruction[i]),"0",1))||(!strncasecmp(&(instruction[i]),"1",1))){
      temp[0]=instruction[i];
      useBeam[count]=atoi(&(temp[0]));
      count--;
    }
  } 
  return;
}/*setBeamsToRead*/


/*####################################################*/
/*read list of beams to skip*/

void setBeamsToSkip(char *useBeam,char *instruction)
{
  int i=0;
  char temp[1];

  /*default is all*/
  for(i=0;i<8;i++)useBeam[i]=1;

  for(i=strlen(instruction)-1;i>=0;i--){
    if((!strncasecmp(&(instruction[i]),"1",1))||(!strncasecmp(&(instruction[i]),"2",1))||\
       (!strncasecmp(&(instruction[i]),"3",1))||(!strncasecmp(&(instruction[i]),"4",1))||\
       (!strncasecmp(&(instruction[i]),"5",1))||(!strncasecmp(&(instruction[i]),"6",1))||\
       (!strncasecmp(&(instruction[i]),"7",1))||(!strncasecmp(&(instruction[i]),"8",1))){
      temp[0]=instruction[i];
      useBeam[atoi(&(temp[0]))-1]=0;
    }
  } 
  return;
}/*setBeamsToSkip*/


/*####################################################*/
/*read list of beams to read*/

void setBeamsToRead(char *useBeam,char *instruction)
{
  int i=0;
  char temp[1];

  /*default is none*/
  for(i=0;i<8;i++)useBeam[i]=0;

  for(i=strlen(instruction)-1;i>=0;i--){
    if((!strncasecmp(&(instruction[i]),"1",1))||(!strncasecmp(&(instruction[i]),"2",1))||\
       (!strncasecmp(&(instruction[i]),"3",1))||(!strncasecmp(&(instruction[i]),"4",1))||\
       (!strncasecmp(&(instruction[i]),"5",1))||(!strncasecmp(&(instruction[i]),"6",1))||\
       (!strncasecmp(&(instruction[i]),"7",1))||(!strncasecmp(&(instruction[i]),"8",1))){
      temp[0]=instruction[i];
      useBeam[atoi(&(temp[0]))-1]=1;
    }
  }
  return;
}/*setBeamsToSkip*/


/*####################################################*/
/*check that number of waves match*/

int checkNwavesDF(int nRead,int nWaves,char *varName)
{
  if(nRead!=nWaves){
    errorf("number of waves mismatch for %s: read %d, expecting %d\n",varName,nRead,nWaves);
    return(-1);
  }

  return(0);
}/*checkNwavesDF*/


/*####################################################*/
/*tidy GEDI HDF data structire*/

gediHDF *tidyGediHDF(gediHDF *hdfData)
{

  if(hdfData){
    TTIDY((void **)hdfData->wave,hdfData->nTypeWaves);
    TTIDY((void **)hdfData->ground,hdfData->nTypeWaves);
    TIDY(hdfData->waveID);
    hdfData->pulse=NULL;     /*as this is repeated in gediIO*/
    TIDY(hdfData->z0);       /*wave top elevations*/
    TIDY(hdfData->zN);       /*wave bottom elevations*/
    TIDY(hdfData->lon);     /*longitudes*/
    TIDY(hdfData->lat);     /*latitudes*/
    TIDY(hdfData->slope);    /*ground slope*/
    TIDY(hdfData->gElev);    /*ground elevation, CofG*/
    TIDY(hdfData->demElev);  /*ground elevation, DEM*/
    TIDY(hdfData->beamDense);/*beam density*/
    TIDY(hdfData->pointDense);/*point density*/
    TIDY(hdfData->zen);      /*scan angles, or mean angles*/
    TIDY(hdfData->sInd);
    TIDY(hdfData->nBins);
    TIDY(hdfData->solarElev);
    TIDY(hdfData);
  }

  return(hdfData);
}/*tidyHDFdata*/


/*###################################################*/
/*read GEDI HDF file*/

dataStruct *unpackHDFgedi(char *namen,gediIOstruct *gediIO,gediHDF **hdfGedi,int numb)
{
  int i=0;
  int sBin=0,eBin=0;
  float zTop=0,maxP=0;
  float *setPulseRange(gediIOstruct *);
  double sOff=0,eOff=0;
  dataStruct *data=NULL;
  void findPCLends(int *,int *,float *,int);


  /*read data from file if needed*/
  if(*hdfGedi==NULL){
    ASSIGN_CHECKNULL_RETNULL(*hdfGedi,readGediHDF(namen,gediIO));
    gediIO->nFiles=hdfGedi[0]->nWaves;
  }/*read data if needed*/

  /*allocate space*/
  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    errorf("error control allocation.\n");
    return(NULL);
  }

  /*copy header*/
  if(hdfGedi[0]->varBins==0)data->nBins=hdfGedi[0]->nBins[0];
  else                      data->nBins=hdfGedi[0]->nBins[numb];
  data->res=fabs(hdfGedi[0]->z0[numb]-hdfGedi[0]->zN[numb])/(float)data->nBins;
  data->nWaveTypes=hdfGedi[0]->nTypeWaves;
  if(data->nWaveTypes<=0)data->nWaveTypes=1;
  data->useType=0;
  data->demGround=0;
  data->pSigma=hdfGedi[0]->pSigma;
  data->fSigma=hdfGedi[0]->fSigma;
  if(hdfGedi[0]->beamDense){
    data->beamDense=hdfGedi[0]->beamDense[numb];
    data->pointDense=hdfGedi[0]->pointDense[numb];
  }
  data->useID=1;
  strcpy(data->waveID,&(hdfGedi[0]->waveID[numb*hdfGedi[0]->idLength]));
  data->lon=hdfGedi[0]->lon[numb];
  data->lat=hdfGedi[0]->lat[numb];
  data->zen=hdfGedi[0]->zen[numb];

  if(gediIO->ground){
    data->slope=hdfGedi[0]->slope[numb];
    data->gElev=hdfGedi[0]->slope[numb];
  }
  data->usable=1;

  /*find start and end of waveform if using PCL*/
  if(gediIO->pcl||gediIO->pclPhoton){
    findPCLends(&sBin,&eBin,&hdfGedi[0]->wave[data->useType][hdfGedi[0]->sInd[numb]],data->nBins);
    /*elevation offsets*/
    sOff=(double)sBin*(double)data->res;
    eOff=(double)(data->nBins-eBin)*(double)data->res;
    /*trim bins*/
    data->nBins=eBin-sBin;
  }else{
    sBin=0;
    eBin=data->nBins;
    sOff=eOff=0.0;
  }

  /*point to arrays rather than copy*/
  ASSIGN_CHECKNULL_RETNULL(data->wave,fFalloc(data->nWaveTypes,"waveform",0));
  ASSIGN_CHECKNULL_RETNULL(data->wave[0],falloc((uint64_t)data->nBins,"waveform",0));
  memcpy(data->wave[0],&(hdfGedi[0]->wave[data->useType][hdfGedi[0]->sInd[numb]+sBin]),data->nBins*4);

  if(gediIO->ground){
    ASSIGN_CHECKNULL_RETNULL(data->ground,fFalloc(data->nWaveTypes,"ground waveform",0));
    ASSIGN_CHECKNULL_RETNULL(data->ground[0],falloc((uint64_t)data->nBins,"waveform",0));
    memcpy(data->ground[0],&(hdfGedi[0]->ground[data->useType][hdfGedi[0]->sInd[numb]+sBin]),data->nBins*4);
  }else data->ground=NULL;

  /*read pulse*/
  if((hdfGedi[0]->nPbins>0)&&(gediIO->pulse==NULL)){
    if(!(gediIO->pulse=(pulseStruct *)calloc(1,sizeof(pulseStruct)))){
      errorf("error pulse allocation.\n");
      return(NULL);
    }
    gediIO->pulse->y=hdfGedi[0]->pulse;
    gediIO->pulse->nBins=hdfGedi[0]->nPbins;
    gediIO->pRes=gediIO->pulse->pRes=hdfGedi[0]->pRes;
    ASSIGN_CHECKNULL_RETNULL(gediIO->pulse->x,setPulseRange(gediIO));

    /*allocate denoising structure if needed*/
    if(gediIO->den==NULL){
      if(!(gediIO->den=(denPar *)calloc(1,sizeof(denPar)))){
        errorf("error denoising parameter allocation.\n");
        return(NULL);
      }
    }

    /*if doing PCL, find peak frequency*/
    if((gediIO->pcl||gediIO->pclPhoton)&&(gediIO->pulse!=NULL)){
      ISINTRETNULL(setPeakChirp(gediIO->pulse));
    }/*peak frequency if needed*/

    /*and copy in to denoising structure*/
    gediIO->den->pBins=gediIO->pulse->nBins;
    gediIO->den->pulse=fFalloc(2,"deconPulse",0);
    ASSIGN_CHECKNULL_RETNULL(gediIO->den->pulse[0],falloc(gediIO->den->pBins,"deconPulse",1));
    ASSIGN_CHECKNULL_RETNULL(gediIO->den->pulse[1],falloc(gediIO->den->pBins,"deconPulse",2));
    memcpy(gediIO->den->pulse[0],gediIO->pulse->x,sizeof(float)*gediIO->den->pBins);
    memcpy(gediIO->den->pulse[1],gediIO->pulse->y,sizeof(float)*gediIO->den->pBins);
    ASSIGN_CHECKNULL_RETNULL(gediIO->den->matchPulse,falloc(gediIO->den->pBins,"matchPulse",0));
    memcpy(gediIO->den->matchPulse,gediIO->pulse->y,sizeof(float)*gediIO->den->pBins);
    maxP=-10000.0;
    for(i=0;i<gediIO->den->pBins;i++){
      if(gediIO->den->pulse[1][i]>maxP){
        maxP=gediIO->den->pulse[1][i];
        gediIO->den->maxPbin=i;
      }
    }
  }else if(hdfGedi[0]->nPbins==0){
    gediIO->pulse=NULL;
  }/*pulse reading*/

  /*count energy*/
  ASSIGN_CHECKNULL_RETNULL(data->totE,falloc((uint64_t)data->nWaveTypes,"totE",0));
  data->totE[data->useType]=0.0;
  for(i=0;i<data->nBins;i++)data->totE[data->useType]+=data->wave[0][i];

  /*elevation needs making and resolution passing to structures*/
  gediIO->res=data->res;
  if(gediIO->gFit)gediIO->gFit->res=data->res;
  if(gediIO->den)gediIO->den->res=data->res;
  ASSIGN_CHECKNULL_RETNULL(data->z,dalloc(data->nBins,"z",0));
  /*which way up are we?*/
  if(hdfGedi[0]->z0[numb]>hdfGedi[0]->zN[numb])zTop=hdfGedi[0]->z0[numb]-sOff;
  else                                         zTop=hdfGedi[0]->zN[numb]-eOff;

  for(i=0;i<data->nBins;i++)data->z[i]=(double)(zTop-(float)i*data->res);

  /*set up number of messages*/
  if(hdfGedi[0]->nWaves>gediIO->nMessages)gediIO->nMessages=(int)(hdfGedi[0]->nWaves/gediIO->nMessages);
  else                                    gediIO->nMessages=1;

  return(data);
}/*unpackHDFgedi*/


/*####################################################*/
/*find ends of signal for PCL*/

void findPCLends(int *sBin,int *eBin,float *wave,int nBins)
{
  int i=0;
  float max=0,thresh=0;

  /*find maximum absolute value*/
  max=-1000.0;
  for(i=0;i<nBins;i++){
    if(fabs(wave[i])>max)max=fabs(wave[i]);
  }
  thresh=max*0.000001;

  /*find the start*/
  *sBin=0;
  for(i=0;i<nBins;i++){
    if(fabs(wave[i])>thresh){
      *sBin=i-1;
      break;
    }
  }
  if(*sBin<0)*sBin=0;

  /*find the end*/
  *eBin=nBins;
  for(i=nBins-1;i>=0;i--){
    if(fabs(wave[i])>thresh){
      *eBin=i+1;
      break;
    } 
  }
  if(*eBin>nBins)*eBin=nBins;

  return;
}/*findPCLends*/


/*####################################################*/
/*set range from pulse file*/

float *setPulseRange(gediIOstruct *gediIO)
{
  int i=0;
  float *x=NULL;
  float max=0;

  /*allocate space*/
  ASSIGN_CHECKNULL_RETNULL(x,falloc(gediIO->pulse->nBins,"pulse range",0));

  /*assign values and check for max*/
  max=-10000.0;
  for(i=0;i<gediIO->pulse->nBins;i++){
    x[i]=(float)i*gediIO->pRes;

    if((gediIO->pcl==0)&&(gediIO->pclPhoton==0)){
      if(gediIO->pulse->y[i]>max){
        max=gediIO->pulse->y[i];
        gediIO->pulse->centBin=i;
      }
    }
  }

  /*to allow for chirps*/
  if(gediIO->pcl||gediIO->pclPhoton)gediIO->pulse->centBin=(int)(gediIO->pulse->nBins/2);

  return(x);
}/*setPulseRange*/


/*####################################################*/
/*read LVIS HDF file*/

dataStruct *unpackHDFlvis(char *namen,lvisHDF **hdfLvis,gediIOstruct *gediIO,int numb)
{
  int i=0,botBin=0;
  int findLvisBottom(float *wave,int nBins);
  dataStruct *data=NULL;
  float *tempPulse=NULL;
  float *pulseLenFromTX(float *,int);
  double dx=0,dy=0,scale=0;  /*for padded LVIS files*/

  /*read data if needed*/
  if(*hdfLvis==NULL){
    ASSIGN_CHECKNULL_RETNULL(*hdfLvis,readLVIShdf(namen));
    gediIO->nFiles=hdfLvis[0]->nWaves;
    gediIO->ground=0;
  }

  /*allocate space*/
  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    errorf("error control allocation.\n");
    return(NULL);
  }
  data->useID=1;
  data->nBins=hdfLvis[0]->nBins;
  data->nWaveTypes=1;
  data->useType=0;
  ASSIGN_CHECKNULL_RETNULL(data->wave,fFalloc(data->nWaveTypes,"waveform",0));
  ASSIGN_CHECKNULL_RETNULL(data->wave[0],falloc((uint64_t)data->nBins,"waveform",0));
  ASSIGN_CHECKNULL_RETNULL(data->totE,falloc((uint64_t)data->nWaveTypes,"totE",0));
  ASSIGN_CHECKNULL_RETNULL(data->z,dalloc(data->nBins,"z",0));
  data->ground=NULL;
  data->demGround=0;
  data->pSigma=-1.0;    /*nonesense pulse length*/
  data->fSigma=-1.0;    /*nonesense footprint width*/
  data->usable=1;

  /*copy data to structure*/
  data->zen=hdfLvis[0]->zen[numb];
  data->res=fabs(hdfLvis[0]->z0[numb]-hdfLvis[0]->z1023[numb])/(float)hdfLvis[0]->nBins;
  if(gediIO->den)gediIO->den->res=data->res;
  if(gediIO->gFit)gediIO->gFit->res=data->res;
  data->totE[data->useType]=0.0;
  for(i=0;i<hdfLvis[0]->nBins;i++){
    data->wave[data->useType][i]=(float)hdfLvis[0]->wave[numb][i];
    data->z[i]=(double)(hdfLvis[0]->z0[numb]-(float)i*data->res);
    data->totE[data->useType]+=data->wave[data->useType][i];
  }
  if(gediIO->den->res<TOL)data->usable=0;
  if(data->totE[data->useType]<=0.0)data->usable=0;
  if(gediIO->ground==0){   /*set to blank*/
    data->cov=-1.0;
    data->gLap=-1.0;
    data->gMinimum=-1.0;
    data->gInfl=-1.0;
  }

  /*for setting coordinate, use range of lowest ground return*/
  ASSIGN_CHECKINT_RETNULL(botBin,findLvisBottom(data->wave[data->useType],data->nBins));

  dx=hdfLvis[0]->lon1023[numb]-hdfLvis[0]->lon0[numb];
  dy=hdfLvis[0]->lat1023[numb]-hdfLvis[0]->lat0[numb];
  scale=(double)botBin/(double)hdfLvis[0]->nBins;
  data->lon=hdfLvis[0]->lon0[numb]+scale*dx;
  data->lat=hdfLvis[0]->lat0[numb]+scale*dy;
  data->lfid=hdfLvis[0]->lfid[numb];
  data->shotN=hdfLvis[0]->shotN[numb];
  snprintf(data->waveID, sizeof(data->waveID), "%d.%d", hdfLvis[0]->lfid[numb], hdfLvis[0]->shotN[numb]);

  /*analyse pulse*/
  if(gediIO->readPsigma){
    ASSIGN_CHECKNULL_RETNULL(tempPulse,falloc((uint64_t)hdfLvis[0]->pBins,"temp pulse",0));
    for(i=0;i<hdfLvis[0]->pBins;i++)tempPulse[i]=(float)hdfLvis[0]->pulse[numb][i];
    float *data_pSigma;
    ASSIGN_CHECKNULL_RETNULL(data_pSigma,pulseLenFromTX(tempPulse,hdfLvis[0]->pBins));
    data->pSigma=*data_pSigma;
    TIDY(tempPulse);
  }else data->pSigma=gediIO->pSigma;

  if(data->fSigma<0.0)data->fSigma=gediIO->fSigma;
  if(data->pSigma<0.0)data->pSigma=gediIO->pSigma;

  /*set up number of messages*/
  if((gediIO->nMessages>0)&&(hdfLvis[0]->nWaves>gediIO->nMessages))gediIO->nMessages=(int)(hdfLvis[0]->nWaves/gediIO->nMessages);
  else                                 gediIO->nMessages=1;

  return(data);
}/*unpackHDFlvis*/


/*####################################################*/
/*find bottom bin for determing coordinate*/

int findLvisBottom(float *wave,int nBins)
{
  int i=0;
  float *tempWave=NULL;
  denPar den;
  char found=0;
  void setDenoiseDefault(denPar *);

  /*set denoising parameters*/
  setDenoiseDefault(&den);
  den.varNoise=1;
  den.statsLen=10.0;
  den.noiseTrack=0;
  den.threshScale=4.0;

  /*find point to calcuklate coordinate from*/
  ASSIGN_CHECKNULL_RETINT(tempWave,processFloWave(wave,nBins,&den,1.0));
  found=0;
  for(i=nBins-1;i>=0;i--){
    if(tempWave[i]>TOL){
      found=1;
      break;
    }
  }
  TIDY(tempWave);
  if(found==0)i=nBins/2;

  return(i);
}/*findLvisBottom*/


/*####################################################*/
/*determine pulse width from TXwave*/

float *pulseLenFromTX(float *pulse,int nBins)
{
  int i=0;
  float *pSigma=0;
  float *denoised=NULL;
  float tot=0,CofG=0;
  denPar den;
  void setDenoiseDefault(denPar *);

  /*denoise*/
  setDenoiseDefault(&den);
  den.varNoise=1;
  den.threshScale=5.0;
  den.noiseTrack=1;
  den.minWidth=5;
  den.statsLen=3.0;
  den.res=0.15;
  ASSIGN_CHECKNULL_RETNULL(denoised,processFloWave(pulse,nBins,&den,1.0));

  /*CofG*/
  CofG=tot=0.0;
  for(i=0;i<nBins;i++){
    CofG+=(float)i*den.res*denoised[i];
    tot+=denoised[i];
  }
  if(tot>0.0){
    CofG/=tot;

    *pSigma=0.0;
    for(i=0;i<nBins;i++)*pSigma+=((float)i*den.res*-CofG)*((float)i*den.res-CofG)*denoised[i];
    *pSigma=sqrt(*pSigma/tot);

  }else *pSigma=-1.0;

  TIDY(denoised);
  return(pSigma);
}/*pulseLenFromTX*/


/*####################################################*/
/*read LVIS binary data*/

dataStruct *readBinaryLVIS(char *namen,lvisLGWstruct *lvis,int numb,gediIOstruct *gediIO)
{
  int i=0,botBin=0;
  int findLvisBottom(float *,int);
  dataStruct *data=NULL;
  float *tempPulse=NULL;
  float *pulseLenFromTX(float *,int);
  double dx=0,dy=0,scale=0;

  /*do we need to read all the data*/
  if(lvis->data==NULL){
    /*read data*/
    ASSIGN_CHECKNULL_RETNULL(lvis->data,readLVISlgw(namen,lvis));
    gediIO->nFiles=lvis->nWaves;
  }


  /*allocate space*/
  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    errorf("error control allocation.\n");
    return(NULL);
  }
  data->useID=1;
  data->nBins=lvis->nBins;
  data->useType=0;
  data->nWaveTypes=1;
  ASSIGN_CHECKNULL_RETNULL(data->wave,fFalloc(data->nWaveTypes,"waveform",0));
  ASSIGN_CHECKNULL_RETNULL(data->wave[data->useType],falloc((uint64_t)data->nBins,"waveform",0));
  ASSIGN_CHECKNULL_RETNULL(data->z,dalloc(data->nBins,"z",0));
  data->ground=NULL;
  data->demGround=0;
  data->pSigma=-1.0;    /*nonesense pulse length*/
  data->fSigma=-1.0;    /*nonesense footprint width*/
  data->usable=1;
  ASSIGN_CHECKNULL_RETNULL(data->totE,falloc((uint64_t)data->nWaveTypes,"totE",0));

  /*copy data to structure*/
  data->zen=lvis->data[numb].zen;
  data->res=(lvis->data[numb].z0-lvis->data[numb].z431)/(float)lvis->nBins;
  if(gediIO->den)gediIO->den->res=data->res;
  if(gediIO->gFit)gediIO->gFit->res=data->res;
  data->totE[data->useType]=0.0;
  for(i=0;i<lvis->nBins;i++){
    data->wave[data->useType][i]=(float)lvis->data[numb].rxwave[i];
    data->z[i]=(double)(lvis->data[numb].z0-(float)i*data->res);
    data->totE[data->useType]+=data->wave[data->useType][i];
  }
  if(gediIO->den->res<TOL)data->usable=0;
  if(data->totE[data->useType]<=0.0)data->usable=0;
  if(gediIO->ground==0){   /*set to blank*/
    data->cov=-1.0;
    data->gLap=-1.0;
    data->gMinimum=-1.0;
    data->gInfl=-1.0;
  }

  /*find point to calculate coordinate from*/
  ASSIGN_CHECKINT_RETNULL(botBin,findLvisBottom(data->wave[data->useType],data->nBins));
  dx=lvis->data[numb].lon431-lvis->data[numb].lon0;
  dy=lvis->data[numb].lat431-lvis->data[numb].lat0;
  scale=(double)botBin/432.0;
  data->lon=lvis->data[numb].lon0+dx*scale;
  data->lat=lvis->data[numb].lat0+dy*scale;
  snprintf(data->waveID, sizeof(data->waveID), "%d.%d", lvis->data[numb].lfid, lvis->data[numb].shotN);

  /*analyse pulse*/
  if(gediIO->readPsigma){
    ASSIGN_CHECKNULL_RETNULL(tempPulse,falloc(80,"temp pulse",0));
    for(i=0;i<80;i++)tempPulse[i]=(float)lvis->data[numb].rxwave[i];
    float *data_pSigma;
    ASSIGN_CHECKNULL_RETNULL(data_pSigma,pulseLenFromTX(tempPulse,80));
    data->pSigma=*data_pSigma;
    TIDY(tempPulse);
  }else data->pSigma=gediIO->pSigma;

  if(data->fSigma<0.0)data->fSigma=gediIO->fSigma;
  if(data->pSigma<0.0)data->pSigma=gediIO->pSigma;


  /*set up number of messages*/
  if((gediIO->nMessages>0)&&(lvis->nWaves>gediIO->nMessages))gediIO->nMessages=(int)(lvis->nWaves/gediIO->nMessages);
  else                              gediIO->nMessages=1;

  return(data);
}/*readBinaryLVIS*/


/*####################################*/
/*read lasFile and save relevant data*/

pCloudStruct *readALSdata(lasFile *las,gediRatStruct *gediRat,int nFile)
{
  int j=0;
  uint32_t i=0;
  uint32_t pUsed=0;    /*number of points used*/
  float decThresh=0;
  double x=0,y=0,z=0;
  pCloudStruct *data=NULL;
  char hasWave=0;   /*has waveform data, to save RAM*/
  char useFile=0,usePoint=0;
  char checkMultiFiles(lasFile *,int,double **,double);
  char checkMultiPoints(double,double,double,int,double **,double);

  /*allocate maximum number of points*/
  if(!(data=(pCloudStruct *)calloc(1,sizeof(pCloudStruct)))){
    errorf("error pCloudStruct allocation.\n");
    return(NULL);
  }

  /*set nonsense bounds*/
  data->bounds[0]=data->bounds[1]=data->bounds[2]=10000000000.0;
  data->bounds[3]=data->bounds[4]=data->bounds[5]=-10000000000.0;

  /*is file needed?*/
  if(!gediRat->readALSonce)useFile=checkFileBounds(las,gediRat->globMinX,gediRat->globMaxX,gediRat->globMinY,gediRat->globMaxY);
  else                     useFile=checkMultiFiles(las,gediRat->gNx,gediRat->coords,gediRat->maxSep);

  /*check file is needed*/
  if(useFile){
    ASSIGN_CHECKNULL_RETNULL(data->x,dalloc(las->nPoints,"x",0));
    ASSIGN_CHECKNULL_RETNULL(data->y,dalloc(las->nPoints,"y",0));
    ASSIGN_CHECKNULL_RETNULL(data->z,dalloc(las->nPoints,"z",0));
    ASSIGN_CHECKNULL_RETNULL(data->refl,ialloc(las->nPoints,"refl",0));
    ASSIGN_CHECKNULL_RETNULL(data->class,uchalloc((uint64_t)las->nPoints,"class",0));
    ASSIGN_CHECKNULL_RETNULL(data->nRet,challoc((uint64_t)las->nPoints,"nRet",0));
    ASSIGN_CHECKNULL_RETNULL(data->retNumb,challoc((uint64_t)las->nPoints,"retNumb",0));
    if(!(data->scanAng=(int16_t *)calloc(las->nPoints,sizeof(int16_t)))){
      errorf("error in input filename structure.\n");
      return(NULL);
    }
    ASSIGN_CHECKNULL_RETNULL(data->packetDes,uchalloc((uint64_t)las->nPoints,"packetDes",0));
    ASSIGN_CHECKNULL_RETNULL(data->grad,fFalloc(las->nPoints,"grad",0));
    for(i=0;i<las->nPoints;i++) {ASSIGN_CHECKNULL_RETNULL(data->grad[i],falloc(3,"grad",i+1));}
    ASSIGN_CHECKNULL_RETNULL(data->time,falloc((uint64_t)las->nPoints,"time",0));
    if(!(data->waveMap=(uint64_t *)calloc(las->nPoints,sizeof(uint64_t)))){
      errorf("error in input filename structure.\n");
      return(NULL);
    }
    if(!(data->waveLen=(uint32_t *)calloc(las->nPoints,sizeof(uint32_t)))){
      errorf("error in input filename structure.\n");
      return(NULL);
    }

    /*loop over points*/
    pUsed=0;
    hasWave=0;
    for(i=0;i<las->nPoints;i++){
      /*read one point*/
      readLasPoint(las,i);
      setCoords(&x,&y,&z,las);

      /*if the point is of use?*/
      if(fabs((float)las->scanAng)<=gediRat->maxScanAng){
        if(!gediRat->readALSonce){
          if((x>=gediRat->globMinX)&&(x<=gediRat->globMaxX)&&(y>=gediRat->globMinY)&&(y<=gediRat->globMaxY)&&\
             (z>-10000.0)&&(z<10000.0))usePoint=1;
          else usePoint=0;
        }else usePoint=checkMultiPoints(x,y,z,gediRat->gNx,gediRat->coords,gediRat->maxSep);
      }else usePoint=0;

      /*are we decimating*/
      if(usePoint){
        if(gediRat->decimate<1.0){
          if(decThresh>gediRat->decimate)usePoint=0;   /*are we accepting this point*/
          if(las->retNumb==las->nRet)decThresh=frand();  /*if last return, draw a new random number*/
        }
      }

      /*if we need to use point*/
      if(usePoint){
        data->x[pUsed]=x;
        data->y[pUsed]=y;
        data->z[pUsed]=z;
        if(las->refl>0)data->refl[pUsed]=(int)las->refl;
        else           data->refl[pUsed]=1;
        data->class[pUsed]=las->classif;
        data->nRet[pUsed]=(char)las->nRet;
        data->retNumb[pUsed]=(char)las->retNumb;
        data->scanAng[pUsed]=las->scanAng;

        /*determine data bounds*/
        if(x<data->bounds[0])data->bounds[0]=x;
        if(y<data->bounds[1])data->bounds[1]=y;
        if(z<data->bounds[2])data->bounds[2]=z;
        if(x>data->bounds[3])data->bounds[3]=x;
        if(y>data->bounds[4])data->bounds[4]=y;
        if(z>data->bounds[5])data->bounds[5]=z;

        /*record waveform if needed*/
        if(checkOneWave(las)){
          hasWave=1;
          data->packetDes[pUsed]=las->packetDes;
          for(j=0;j<3;j++)data->grad[pUsed][j]=las->grad[j];
          data->time[pUsed]=las->time;
          data->waveMap[pUsed]=las->waveMap;
          data->waveLen[pUsed]=las->waveLen;
        }else{
          data->packetDes[pUsed]=0;
          data->grad[pUsed][0]=data->grad[pUsed][1]=data->grad[pUsed][2]=0.0;
        }

        /*map to octree if needed*/
        if(gediRat->useOctree) {ISINTRETNULL(fillOctree(x,y,z,nFile,pUsed,gediRat->octree));}

        /*count points here*/
        pUsed++;
      }
    }/*point loop*/

    /*trim data arrays*/
    data->nPoints=pUsed;
    if(pUsed>0){
      if(!(data->x=(double *)realloc(data->x,data->nPoints*sizeof(double)))){
        errorf("Balls\n");
        return(NULL);
      }
      if(!(data->y=(double *)realloc(data->y,data->nPoints*sizeof(double)))){
        errorf("Balls\n");
        return(NULL);
      }
      if(!(data->z=(double *)realloc(data->z,data->nPoints*sizeof(double)))){
        errorf("Balls\n");
        return(NULL);
      }
      if(!(data->refl=(int *)realloc(data->refl,data->nPoints*sizeof(int)))){
        errorf("Balls\n");
        return(NULL);
      }
      if(!(data->class=(unsigned char *)realloc(data->class,data->nPoints*sizeof(unsigned char)))){
        errorf("Balls\n");
        return(NULL);
      }
      if(!(data->nRet=(char *)realloc(data->nRet,data->nPoints*sizeof(char)))){
        errorf("Balls\n");
        return(NULL);
      }
      if(!(data->retNumb=(char *)realloc(data->retNumb,data->nPoints*sizeof(char)))){
        errorf("Balls\n");
        return(NULL);
      }
      if(!(data->scanAng=(int16_t *)realloc(data->scanAng,data->nPoints*sizeof(int16_t)))){
        errorf("Balls\n");
        return(NULL);
      }
      if(gediRat->useShadow){
        for(i=data->nPoints;i<las->nPoints-data->nPoints;i++)TIDY(data->grad[i]);
        if(!(data->grad=(float **)realloc(data->grad,data->nPoints*sizeof(float *)))){
          errorf("Balls\n");
          return(NULL);
        }
      }else if(hasWave==0){
        TTIDY((void **)data->grad,las->nPoints);
        data->grad=NULL;
      }
    }else{
      TIDY(data->x);
      TIDY(data->y);
      TIDY(data->z);
      TIDY(data->refl);
      TIDY(data->class);
      TTIDY((void **)data->grad,las->nPoints);
      data->grad=NULL;
    }
    if(hasWave==1){
      data->waveStart=las->waveStart;
      if(!(data->packetDes=(unsigned char *)realloc(data->packetDes,data->nPoints*sizeof(unsigned char)))){
        errorf("Balls\n");
        return(NULL);
      }
      if(!(data->time=(float *)realloc(data->time,data->nPoints*sizeof(float)))){
        errorf("Balls\n");
        return(NULL);
      }
      if(!(data->waveMap=(uint64_t *)realloc(data->waveMap,data->nPoints*sizeof(uint64_t)))){
        errorf("Balls\n");
        return(NULL);
      }
      if(!(data->waveLen=(uint32_t *)realloc(data->waveLen,data->nPoints*sizeof(uint32_t)))){
        errorf("Balls\n");
        return(NULL);
      }
      for(i=data->nPoints;i<las->nPoints-data->nPoints;i++)TIDY(data->grad[i]);
      if(!(data->grad=(float **)realloc(data->grad,data->nPoints*sizeof(float *)))){
        errorf("Balls\n");
        return(NULL);
      }
    }else{  /*clear out all the waveform bits*/
      TIDY(data->packetDes);
      if(!gediRat->useShadow){
        TTIDY((void **)data->grad,las->nPoints);
        data->grad=NULL;
      }
      TIDY(data->time);
      TIDY(data->waveMap);
      TIDY(data->waveLen);
    }
  }else{/*file bounds check*/
    data->nPoints=0;
    data->nRet=NULL;
    data->retNumb=NULL;
    data->packetDes=NULL;
    data->grad=NULL;
    data->time=NULL;
    data->waveMap=NULL;
    data->waveLen=NULL;
    data->x=NULL;
    data->y=NULL;
    data->z=NULL;
    data->refl=NULL;
  }

  data->hasWave=hasWave;
  if(gediRat->readWave&&hasWave){  /*only leave files open for waveform*/
    data->ipoo=las->ipoo;
    las->ipoo=NULL;
  }else{
    data->ipoo=NULL;
  }

  return(data);
}/*readALSdata*/


/*##########################################################################*/
/*see if we need to use this file when batch processing ALS data*/

char checkMultiFiles(lasFile *las,int nCoords,double **coords,double maxSep)
{
  int i=0;
  double maxX=0,minX=0;
  double maxY=0,minY=0;
  char useFile=0;

  useFile=0;
  for(i=0;i<nCoords;i++){
    minX=coords[i][0]-maxSep;
    maxX=coords[i][0]+maxSep;
    minY=coords[i][1]-maxSep;
    maxY=coords[i][1]+maxSep;
    if((las->minB[0]<=maxX)&&(las->minB[1]<=maxY)&&(las->maxB[0]>=minX)&&(las->maxB[1]>=minY)){
      useFile=1;
      break;
    }
  }/*footprint loop*/

  return(useFile);
}/*checkMultiFiles*/


/*##########################################################################*/
/*see if we need to use this point when batch processing ALS data*/

char checkMultiPoints(double x,double y,double z,int nCoords,double **coords,double maxSep)
{
  int i=0;
  double maxX=0,minX=0;
  double maxY=0,minY=0;
  char usePoint=0;

  usePoint=0;
  for(i=0;i<nCoords;i++){
    minX=coords[i][0]-maxSep;
    maxX=coords[i][0]+maxSep;
    minY=coords[i][1]-maxSep;
    maxY=coords[i][1]+maxSep;

    if((x>=minX)&&(x<=maxX)&&(y>=minY)&&(y<=maxY)){
      usePoint=1;
      break;
    }
  }/*footprint loop*/

  return(usePoint);
}/*checkMultiPoints*/


/*####################################*/
/*set GEDI grid or batch*/

int setGediGrid(gediIOstruct *gediIO,gediRatStruct *gediRat)
{
  int readFeetList(gediRatStruct *);
  void setRatBounds(gediRatStruct *);
  int readWavefront(gediRatStruct *,gediIOstruct *);

  /*footprint width*/
  if(gediRat->defWfront==0){   /*regular footprint*/
    if(gediRat->topHat==0)gediRat->maxSep=determineGaussSep(gediIO->fSigma,gediRat->iThresh);
    else                  gediRat->maxSep=gediIO->fSigma;
  }else{    /*read assymetric footprint*/
    ISINTRETINT(readWavefront(gediRat,gediIO));
  }/*footprint width setting*/

  if(gediRat->doGrid){  /*it is a grid*/
    /*number of footprints*/
    gediRat->gNx=(int)((gediRat->gMaxX-gediRat->gMinX)/gediRat->gRes+1);
    gediRat->gNy=(int)((gediRat->gMaxY-gediRat->gMinY)/gediRat->gRes+1);

    /*global bounds*/
    gediRat->globMinX=gediRat->gMinX-gediRat->maxSep;
    gediRat->globMaxX=gediRat->gMaxX+gediRat->maxSep;
    gediRat->globMinY=gediRat->gMinY-gediRat->maxSep;
    gediRat->globMaxY=gediRat->gMaxY+gediRat->maxSep;
  }else if(gediRat->readALSonce){ /*it is a batch*/
    /*read list of coords*/
    if(gediRat->coords==NULL)
      ISINTRETINT(readFeetList(gediRat));
    setRatBounds(gediRat);
  }else{   /*single footprint*/
    gediRat->gNx=gediRat->gNy=1;
    gediRat->globMinX=gediRat->coord[0]-gediRat->maxSep;
    gediRat->globMaxX=gediRat->coord[0]+gediRat->maxSep;
    gediRat->globMinY=gediRat->coord[1]-gediRat->maxSep;
    gediRat->globMaxY=gediRat->coord[1]+gediRat->maxSep;
    ASSIGN_CHECKNULL_RETINT(gediRat->waveIDlist,chChalloc(1,"wave ID list",0));
    ASSIGN_CHECKNULL_RETINT(gediRat->waveIDlist[0],challoc(300,"wave ID list",1));
    snprintf(gediRat->waveIDlist[0], 300, "BEAM.x.%.2f.y.%.2f", gediRat->coord[0], gediRat->coord[1]);
  }

  if((gediIO->nMessages>1)&&(gediRat->gNx*gediRat->gNy)>gediIO->nMessages)gediIO->nMessages=(int)(gediRat->gNx*gediRat->gNy/gediIO->nMessages);
  else                                             gediIO->nMessages=1;

  /*allocate octree if needed*/
  if(gediRat->useOctree){
    ASSIGN_CHECKNULL_RETINT(gediRat->octree,allocateOctree(gediRat->octLevels,gediRat->nOctTop,\
            gediRat->globMinX,gediRat->globMaxX,gediRat->globMinY,gediRat->globMaxY));
  }else gediRat->octree=NULL;

  return(0);
}/*setGediGrid*/


/*###################################################*/

int readWavefront(gediRatStruct *gediRat,gediIOstruct *gediIO)
{
  int i=0,j=0,maxI=0;
  float len=0; /*total=0*/
  char line[20000];
  char *token=NULL;
  void setWavefrontRes(wFrontStruct *,float);
  FILE *ipoo=NULL;

  /*open file*/
  if((ipoo=fopen(gediRat->wavefront->frontFile,"r"))==NULL){
    errorf("Error opening wavefront file \"%s\"\n",gediRat->wavefront->frontFile);
    return(-1);
  }

  /*find file size*/
  j=0;
  maxI=0;
  while(fgets(line,20000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      i=0;
      token=strtok(line,",");
      while(token){
        token=strtok(NULL,",");
        i++;
      }
      if(i>maxI)maxI=i;
      j++;
    }
  }/*file size counting*/
  gediRat->wavefront->nX=maxI;
  gediRat->wavefront->nY=j;

  /*allocate space*/
  ASSIGN_CHECKNULL_RETINT(gediRat->wavefront->front,fFalloc(gediRat->wavefront->nX,"wavefront",0));
  for(i=0;i<gediRat->wavefront->nX;i++){
    ASSIGN_CHECKNULL_RETINT(gediRat->wavefront->front[i],falloc((uint64_t)gediRat->wavefront->nY,"wavefront",j));
    for(j=0;j<gediRat->wavefront->nY;j++)gediRat->wavefront->front[i][j]=-1.0;  /*mark as blank*/
  }/*allocation*/

  /*rewind*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    errorf("fseek error\n");
    return(-1);
  }

  /*read data*/
  /*total=0.0;*/
  j=gediRat->wavefront->nY-1;
  while(fgets(line,20000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      i=0;
      token=strtok(line,",");
      while(token){
        if(strncasecmp(token,",",1)){  /*is there a result*/
          gediRat->wavefront->front[i][j]=atof(token);
          /*total+=gediRat->wavefront->front[i][j];*/
        }
        token=strtok(NULL,",");
        i++;
      }
      j--;
    }
  }/*data reading*/

  /*tidy up*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  /*normalise*/
  /*for(i=0;i<gediRat->wavefront->nX;i++){
    for(j=0;j<gediRat->wavefront->nY;j++){
      if(gediRat->wavefront->front[i][j]>0.0)gediRat->wavefront->front[i][j]/=total;
    }
  }*//*normalisation*/

  /*determine resolution from footprint width*/
  setWavefrontRes(gediRat->wavefront,gediIO->fSigma);

  /*set bounds for reading*/
  gediRat->maxSep=-100000.0;
  len=(gediRat->wavefront->x0>(gediRat->wavefront->nX/2))?(float)gediRat->wavefront->x0*gediRat->wavefront->res:\
                  (float)(gediRat->wavefront->nX-gediRat->wavefront->x0)*gediRat->wavefront->res;
  gediRat->maxSep=(len>gediRat->maxSep)?len:gediRat->maxSep;
  len=(gediRat->wavefront->y0>(gediRat->wavefront->nY/2))?(float)gediRat->wavefront->y0*gediRat->wavefront->res:\
                  (float)(gediRat->wavefront->nY-gediRat->wavefront->y0)*gediRat->wavefront->res;
  gediRat->maxSep=(len>gediRat->maxSep)?len:gediRat->maxSep;

  return(0);
}/*readWavefront*/


/*##########################################################*/
/*determine wavefront resolution and peak*/

void setWavefrontRes(wFrontStruct *wavefront,float fSigma)
{
  int i=0,j=0;
  int ii=0,jj=0;
  int window=0,contN=0;
  float max=0,tot=0,mean=0;
  float xStdev=0,yStdev=0;

  /*find centre*/
  window=2;
  max=-1.0;
  for(i=0;i<wavefront->nX;i++){
    for(j=0;j<wavefront->nY;j++){
      mean=0.0;
      contN=0;
      for(ii=i-window;ii<(i+window);ii++){
        if((ii<0)||(ii>=wavefront->nX))continue;
        for(jj=j-window;jj<(j+window);jj++){
          if((jj<0)||(jj>=wavefront->nY))continue;
          if(wavefront->front[i][j]>=0.0){
            mean+=wavefront->front[ii][jj];
            contN++;
          }
        }
      }
      if(contN>0){
        mean/=(float)contN;
        if(mean>max){
          max=mean;
          wavefront->x0=i;
          wavefront->y0=j;
        }
      }
    }
  }

  /*determine width in two axes*/
  tot=xStdev=0.0;
  for(i=0;i<wavefront->nX;i++){
    if(wavefront->front[i][wavefront->y0]>=0.0){
      xStdev+=pow((float)(i-wavefront->x0),2.0)*wavefront->front[i][wavefront->y0];
      tot+=wavefront->front[i][wavefront->y0];
    }
  }
  xStdev=sqrt(xStdev/tot);
  tot=yStdev=0.0;
  for(j=0;j<wavefront->nY;j++){
    if(wavefront->front[wavefront->x0][j]>=0.0){
      yStdev+=pow((float)(j-wavefront->y0),2.0)*wavefront->front[wavefront->x0][j];
      tot+=wavefront->front[wavefront->x0][j];
    }
  }
  yStdev=sqrt(yStdev/tot);
  wavefront->res=2.0*fSigma/(xStdev+yStdev);

  return;
}/*setWavefrontRes*/


/*####################################*/
/*set bounds from list of coords*/

void setRatBounds(gediRatStruct *gediRat)
{
  int i=0;
  double minX=0,maxX=0;
  double minY=0,maxY=0;

  minX=minY=1000000000.0;
  maxX=maxY=-1000000000.0;
  for(i=0;i<gediRat->gNx;i++){
    if(gediRat->coords[i][0]<minX)minX=gediRat->coords[i][0];
    if(gediRat->coords[i][1]<minY)minY=gediRat->coords[i][1];
    if(gediRat->coords[i][0]>maxX)maxX=gediRat->coords[i][0];
    if(gediRat->coords[i][1]>maxY)maxY=gediRat->coords[i][1];
  }

  gediRat->globMinX=minX-gediRat->maxSep;
  gediRat->globMaxX=maxX+gediRat->maxSep;
  gediRat->globMinY=minY-gediRat->maxSep;
  gediRat->globMaxY=maxY+gediRat->maxSep;

  return;
}/*setRatBounds*/


/*####################################*/
/*read list of coordinates*/

int readFeetList(gediRatStruct *gediRat)
{
  int i=0;
  char line[200],temp1[50];
  char temp2[50],temp3[100];
  char temp4[50],temp5[50];
  FILE *ipoo=NULL;

  /*open file*/
  if((ipoo=fopen(gediRat->coordList,"r"))==NULL){
    errorf("Error opening input file list \"%s\"\n",gediRat->coordList);
    return(-1);
  }


  /*count number of lines*/
  i=0;
  while(fgets(line,200,ipoo)!=NULL)if(strncasecmp(line,"#",1))i++;

  /*allocate space*/
  gediRat->gNx=i;
  gediRat->gNy=1;
  ASSIGN_CHECKNULL_RETINT(gediRat->coords,dDalloc(gediRat->gNx,"coord list",0));
  ASSIGN_CHECKNULL_RETINT(gediRat->waveIDlist,chChalloc(gediRat->gNx,"wave ID list",0));
  ASSIGN_CHECKNULL_RETINT(gediRat->geoCoords,dDalloc(gediRat->gNx,"geolocated coord list",0));

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    errorf("fseek error\n");
    return(-1);
  }

  /*read coordinate list*/
  i=0;
  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      ASSIGN_CHECKNULL_RETINT(gediRat->coords[i],dalloc(2,"coord list",i+1));
      ASSIGN_CHECKNULL_RETINT(gediRat->geoCoords[i],dalloc(2,"geolocated coord list",i+1));
      if(sscanf(line,"%s %s %s %s %s",temp1,temp2,temp3,temp4,temp5)==5){ /*read coord, waveID and shofted coord*/
        gediRat->coords[i][0]=atof(temp1);
        gediRat->coords[i][1]=atof(temp2);
        ASSIGN_CHECKNULL_RETINT(gediRat->waveIDlist[i],challoc((int)strlen(temp3)+1,"wave ID list",i+1));
        strcpy(gediRat->waveIDlist[i],temp3);
        gediRat->geoCoords[i][0]=atof(temp4);
        gediRat->geoCoords[i][1]=atof(temp5);
      }else if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){ /*read coord and waveID*/
        gediRat->coords[i][0]=gediRat->geoCoords[i][0]=atof(temp1);
        gediRat->coords[i][1]=gediRat->geoCoords[i][1]=atof(temp2);
        ASSIGN_CHECKNULL_RETINT(gediRat->waveIDlist[i],challoc((int)strlen(temp3)+1,"wave ID list",i+1));
        strcpy(gediRat->waveIDlist[i],temp3);
      }else if(sscanf(line,"%s %s",temp1,temp2)==2){
        gediRat->coords[i][0]=gediRat->geoCoords[i][0]=atof(temp1);
        gediRat->coords[i][1]=gediRat->geoCoords[i][1]=atof(temp2);
        snprintf(temp3, sizeof(temp3), "%f.%f", gediRat->coords[i][0], gediRat->coords[i][1]);
        ASSIGN_CHECKNULL_RETINT(gediRat->waveIDlist[i],challoc((int)strlen(temp3)+1,"wave ID list",i+1));
        strcpy(gediRat->waveIDlist[i],temp3);
      }else{
        errorf("coord list reading error \"%s\"\n",line);
        return(-1);
      }
      i++;
    }
  }

  /*close file*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(0);
}/*readFeetList*/


/*####################################*/
/*set GEDI pulse*/

int setGediPulse(gediIOstruct *gediIO,gediRatStruct *gediRat)
{
  int i=0;
  float fwhm=0;   /*FWHM in metres*/
  float x=0,y=0;
  float max=0,tot=0;
  int readSimPulse(gediIOstruct *);

  if(!(gediIO->pulse=(pulseStruct *)calloc(1,sizeof(pulseStruct)))){
    errorf("error pulseStruct allocation.\n");
    return(-1);
  }


  if(gediIO->readPulse==0){  /*Gaussian pulse*/
    /*pulse length*/
    /*calculate sigma from FWHM*/
    if(gediIO->pSigma<0.0){  /*GEDI unless specificed*/
      fwhm=gediIO->pFWHM*0.2998/2.0;  /*time for two way*/
      gediIO->pSigma=fwhm/2.35482;  /* =2*sqrt(2*ln2) */
    }

    if(gediIO->pSigma>0.0){  /*if we are using a pulse width*/
      /*determine number of bins*/
      gediIO->pulse->nBins=0;
      x=0.0;
      if(gediRat->iThresh<=0.0)gediRat->iThresh=0.0006;  /*prevent infinite tolerances*/
      do{
        y=(float)gaussian((double)x,(double)gediIO->pSigma,0.0);
        x+=gediIO->pRes;
        gediIO->pulse->nBins+=2;  /*both sides of peak*/
      }while(y>=gediRat->iThresh);

      ASSIGN_CHECKNULL_RETINT(gediIO->pulse->x,falloc((uint64_t)gediIO->pulse->nBins,"pulse x",0));
      ASSIGN_CHECKNULL_RETINT(gediIO->pulse->y,falloc((uint64_t)gediIO->pulse->nBins,"pulse y",0));
      gediIO->pulse->centBin=(int)(gediIO->pulse->nBins/2);

      max=-100.0;
      tot=0.0;
      x=-1.0*(float)gediIO->pulse->centBin*gediIO->pRes;
      for(i=0;i<gediIO->pulse->nBins;i++){
        gediIO->pulse->x[i]=x;
        gediIO->pulse->y[i]=(float)gaussian((double)x,(float)gediIO->pSigma,0.0);
        if(gediIO->pulse->y[i]>max){
          max=gediIO->pulse->y[i];
          gediIO->pulse->centBin=i;
        }
        tot+=gediIO->pulse->y[i];
        x+=gediIO->pRes;
      }
      /*normalise to cope with rounding*/
      for(i=0;i<gediIO->pulse->nBins;i++){
        gediIO->pulse->y[i]/=tot;
      }
    }else{  /*dirac-delta*/
      gediIO->pulse->nBins=1;
      ASSIGN_CHECKNULL_RETINT(gediIO->pulse->x,falloc((uint64_t)gediIO->pulse->nBins,"pulse x",0));
      ASSIGN_CHECKNULL_RETINT(gediIO->pulse->y,falloc((uint64_t)gediIO->pulse->nBins,"pulse y",0));
      gediIO->pulse->centBin=0;

      gediIO->pulse->x[0]=0.0;
      gediIO->pulse->y[0]=1.0;
    }
  }else{  /*read the pulse from a file*/
    ISINTRETINT(readSimPulse(gediIO));
  }
  
  /*determine peak frequency if needed*/
  if((gediIO->pcl||gediIO->pclPhoton)&&(gediIO->pulse!=NULL)){
    setPeakChirp(gediIO->pulse);
  }/*peak frequency if needed*/


  return(0);
}/*setGediPulse*/


/*####################################*/
/*set peak frequency of chirp*/

int setPeakChirp(pulseStruct *pulse)
{
  int i=0,numb=0;
  // int startInd=0;
  //int status=0;
  float maxAmp=0,thresh=0.0;
  float *absData=NULL;
  double *data=NULL;
  char hasStarted=0;


  /*copy data to correct type*/
  numb=pow(2.0,(float)((int)(log((float)pulse->nBins)/log(2.0)+0.5)+1));
  ASSIGN_CHECKNULL_RETINT(data,dalloc(numb*2,"double pulse",0));
  for(i=0;i<pulse->nBins;i++){
    data[2*i]=(double)pulse->y[i];
    data[2*i+1]=0.0;
  }
  for(;i<2*numb;i++)data[i]=0.0;

  /*apply FFT to pulse*/
  gsl_fft_complex_radix2_forward((gsl_complex_packed_array)data,1,(size_t)numb);

  /*get absolute value*/

  ASSIGN_CHECKNULL_RETINT(absData,falloc(numb,"absolute FFT pulse",0));
  maxAmp=0.0;
  for(i=0;i<numb;i++){
    absData[i]=sqrt(data[2*i]*data[2*i]+data[2*i+1]*data[2*i+1]);
    if(absData[i]>maxAmp)maxAmp=absData[i];
  }
  TIDY(data);

  /*determine bin of crossing point*/
  thresh=maxAmp/2.0;
  hasStarted=0;
  for(i=0;i<numb/2;i++){
    if(!hasStarted){
      if(absData[i]>=thresh)hasStarted=1;
      // startInd=i;
    }else{
      if(absData[i]<=thresh)break;
    }
  }
  TIDY(absData);

  /*scale to frequency*/
  pulse->peakFreq=(float)i/(float)(numb/2)*(2.998*pow(10,8))/pulse->pRes;
  msgf("Pulse peak frequency is %.2f MHz\n",pulse->peakFreq/1000000.0);

  return 0;
}/*setPeakChirp*/


/*####################################*/
/*read pulse to use for simulator*/

int readSimPulse(gediIOstruct *gediIO)
{
  int i=0;
  float CofG=0,tot=0,centre=0;
  float minSep=0,max=0;
  char line[400];
  char temp1[100],temp2[100];
  FILE *ipoo=NULL;

  if((ipoo=fopen(gediIO->pulseFile,"r"))==NULL){
    errorf("Error opening input file %s\n",gediIO->pulseFile);
    return(-1);
  }

  /*count number of bins*/
  gediIO->pulse->nBins=0;
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))gediIO->pulse->nBins++;

  ASSIGN_CHECKNULL_RETINT(gediIO->pulse->x,falloc((uint64_t)gediIO->pulse->nBins,"pulse x",0));
  ASSIGN_CHECKNULL_RETINT(gediIO->pulse->y,falloc((uint64_t)gediIO->pulse->nBins,"pulse y",0));

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    errorf("fseek error\n");
    return(-1);
  }

  /*read data*/
  i=0;
  while(fgets(line,400,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s",temp1,temp2)==2){
        gediIO->pulse->x[i]=atof(temp1);
        gediIO->pulse->y[i]=atof(temp2);
        i++;
      }
    }
  }
  gediIO->pRes=fabs(gediIO->pulse->x[gediIO->pulse->nBins-1]-gediIO->pulse->x[0])/(float)(gediIO->pulse->nBins-1);
  gediIO->pulse->pRes=gediIO->pRes;

  /*determine maximum to centre and total to normalise*/
  tot=0.0;
  CofG=0.0;
  max=-1000.0;
  for(i=0;i<gediIO->pulse->nBins;i++){
    CofG+=gediIO->pulse->x[i]*gediIO->pulse->y[i];
    if(gediIO->pulse->y[i]>=max){
      max=gediIO->pulse->y[i];
      centre=gediIO->pulse->x[i];
    }
    tot+=gediIO->pulse->y[i];
  }

  if(tot>0.0)CofG/=tot;
  CofG-=centre;

  if((gediIO->pcl==1)||(gediIO->pclPhoton))centre=gediIO->pulse->x[(int)(gediIO->pulse->nBins/2)];

  /*align pulse*/
  if((gediIO->pcl==0)&&(gediIO->pclPhoton==0)){
    minSep=1000.0;
    gediIO->pSigma=0.0;
    for(i=0;i<gediIO->pulse->nBins;i++){
      gediIO->pulse->x[i]-=centre;

      if(fabs(gediIO->pulse->x[i])<minSep){
        minSep=fabs(gediIO->pulse->x[i]);
        gediIO->pulse->centBin=i;
      }
    }
  }else gediIO->pulse->centBin=(int)(gediIO->pulse->nBins/2);  /*if we are using pulse compressed lidar*/

  /*pulse width*/
  gediIO->pSigma=0.0;
  for(i=0;i<gediIO->pulse->nBins;i++){
    gediIO->pSigma+=(gediIO->pulse->x[i]-CofG)*(gediIO->pulse->x[i]-CofG)*gediIO->pulse->y[i];
  }
  gediIO->pSigma=sqrt(gediIO->pSigma/tot);
  gediIO->linkPsig=gediIO->pSigma;

  /*now normalise*/
  //tot=tot-min*(float)gediIO->pulse->nBins;
  //for(i=0;i<gediIO->pulse->nBins;i++)gediIO->pulse->y[i]=(gediIO->pulse->y[i]-min)/tot;
  for(i=0;i<gediIO->pulse->nBins;i++)gediIO->pulse->y[i]/=tot;

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(0);
}/*readSimPulse*/


/*####################################*/
/*set GEDI footprint*/

int setGediFootprint(gediRatStruct *gediRat,gediIOstruct *gediIO)
{
  int i=0;
  float totE=0;
  float az=0;
  double tX=0,tY=0;
  int intersectOctree(gediRatStruct *);


  /*footprint width*/
  az=gediRat->lobeAng*M_PI/180.0;  /*convert anlge to radians*/

  /*number of lobes and allocate*/
  if(gediRat->sideLobe==0)gediRat->nLobes=1;
  else                    gediRat->nLobes=7;
  if(!(gediRat->lobe=(lobeStruct *)calloc(gediRat->nLobes,sizeof(lobeStruct)))){
    errorf("error lobeStruct allocation.\n");
    return(-1);
  }

  /*central footprint*/
  i=0;
  gediRat->lobe[i].coord[0]=gediRat->coord[0];
  gediRat->lobe[i].coord[1]=gediRat->coord[1];
  gediRat->lobe[i].fSigma=gediIO->fSigma;
  gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);
  if(gediRat->sideLobe==0)gediRat->lobe[i].E=1.0;
  else{  /*include side lobes*/
    totE=1.0+0.0599+0.0731+0.0317+0.0319+0.0167+0.0163;
    i=0;
    gediRat->lobe[i].E=1.0/totE;

    /*first southern lobe*/
    i=1;
    gediRat->lobe[i].E=0.0731/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]-20.0*sin(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]-20.0*cos(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);

    /*first nothern lobe*/
    i=2;
    gediRat->lobe[i].E=0.0599/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]+20.0*sin(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]+20.0*cos(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);

    /*western lobe*/
    i=3;
    gediRat->lobe[i].E=0.0319/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]-20.0*cos(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]-20.0*sin(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);

    /*eastern lobe*/
    i=4;
    gediRat->lobe[i].E=0.0317/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]+20.0*cos(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]+20.0*sin(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);

    /*second southern lobe*/
    i=5;
    gediRat->lobe[i].E=0.0167/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]-30.0*sin(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]-30.0*cos(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);

    /*second northern lobe*/
    i=6;
    gediRat->lobe[i].E=0.0163/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]+30.0*cos(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]+30.0*sin(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);
  }/*side lobe test*/


  /*determine min and max bounds*/
  gediRat->minX=gediRat->minY=100000000000.0;
  gediRat->maxX=gediRat->maxY=-1000000000000.0;
  for(i=0;i<gediRat->nLobes;i++){
    tX=gediRat->lobe[i].coord[0]-sqrt(gediRat->lobe[i].maxSepSq);
    if(tX<gediRat->minX)gediRat->minX=tX;
    tX=gediRat->lobe[i].coord[0]+sqrt(gediRat->lobe[i].maxSepSq);
    if(tX>gediRat->maxX)gediRat->maxX=tX;
    tY=gediRat->lobe[i].coord[1]-sqrt(gediRat->lobe[i].maxSepSq);
    if(tY<gediRat->minY)gediRat->minY=tY;
    tY=gediRat->lobe[i].coord[1]+sqrt(gediRat->lobe[i].maxSepSq);
    if(tY>gediRat->maxY)gediRat->maxY=tY;
  }

  /*grid for normalising sampling*/
  if(gediRat->normCover||gediRat->checkCover){
    gediRat->gridRes=1.5;
    gediRat->gX=(int)((float)(gediRat->maxX-gediRat->minX)/gediRat->gridRes)+2;
    gediRat->gY=(int)((float)(gediRat->maxY-gediRat->minY)/gediRat->gridRes)+2;
    gediRat->g0[0]=gediRat->minX+(double)gediRat->gridRes;
    gediRat->g0[1]=gediRat->minY+(double)gediRat->gridRes;
    ASSIGN_CHECKNULL_RETINT(gediRat->nGrid,ialloc(gediRat->gX*gediRat->gY,"nGrid",0));
    for(i=gediRat->gX*gediRat->gY-1;i>=0;i--)gediRat->nGrid[i]=0;
  }

  /*radius to calculate density within*/
  if(gediRat->topHat==0)gediRat->denseRadSq=gediIO->fSigma*gediIO->fSigma*4.0;
  else                  gediRat->denseRadSq=gediIO->fSigma;
  gediRat->pointDense=gediRat->beamDense=0.0;

  /*determine which octree cells are intresected*/
  if(gediRat->useOctree) {ISINTRETINT(intersectOctree(gediRat));}

  return(0);
}/*setGediFootprint*/


/*#####################################################*/
/*determine which top level octree pixels intersect*/

int intersectOctree(gediRatStruct *gediRat)
{
  int i=0,j=0;
  int minI=0,maxI=0;
  int minJ=0,maxJ=0;
  int *markInt(int,int *,int);

  /*reset counters*/
  gediRat->nOct=0;
  TIDY(gediRat->octList);
  gediRat->octList=NULL;

  /*determine bounds to search*/
  minI=(int)((gediRat->minX-gediRat->octree->minX)/(double)gediRat->octree->res);
  maxI=(int)((gediRat->maxX-gediRat->octree->minX)/(double)gediRat->octree->res+0.5);
  minJ=(int)((gediRat->minY-gediRat->octree->minY)/(double)gediRat->octree->res);
  maxJ=(int)((gediRat->maxY-gediRat->octree->minY)/(double)gediRat->octree->res+0.5);

  /*loop over and test*/
  for(i=minI;i<=maxI;i++){
    if((i<0)||(i>=gediRat->octree->nX))continue;
    for(j=minJ;j<=maxJ;j++){
      if((j<0)||(j>=gediRat->octree->nY))continue;
      ASSIGN_CHECKNULL_RETINT(gediRat->octList,markInt(gediRat->nOct,gediRat->octList,i+j*gediRat->octree->nX));
      gediRat->nOct++;
    }
  }
  return(0);
}/*intersectOctree*/


/*##############################################*/
/*copy waveform into HDF structure*/

int packGEDIhdf(waveStruct *waves,gediHDF *hdfData,int waveNumb,gediIOstruct *gediIO,gediRatStruct *gediRat,int *hdfCount,char useID,char *waveID)
{
  int i=0,j=0,start=0,numb=0;
  int nBins=0,idLength=0;
  float *tot=NULL,*cumul=NULL;
  float *thresh=NULL,buff=0;
  char thisWaveID[300];

  numb=*hdfCount;


  /*trim waveform*/
  if(gediIO->pcl==0){
    buff=30.0;
    /*find energies*/
    ASSIGN_CHECKNULL_RETINT(tot,falloc((uint64_t)hdfData->nTypeWaves,"tot",0));
    ASSIGN_CHECKNULL_RETINT(cumul,falloc((uint64_t)hdfData->nTypeWaves,"cumul",0));
    for(j=0;j<hdfData->nTypeWaves;j++){
      tot[j]=cumul[j]=0.0;
      for(i=0;i<waves->nBins;i++)tot[j]+=waves->wave[j][i];
    }


    /*set threshold*/
    ASSIGN_CHECKNULL_RETINT(thresh,falloc((uint64_t)hdfData->nTypeWaves,"thresh",0));
    for(j=0;j<hdfData->nTypeWaves;j++)thresh[j]=0.01*tot[j];
    TIDY(tot);

    /*find waveform start*/
    start=-1;
    for(i=0;i<waves->nBins;i++){
      for(j=0;j<hdfData->nTypeWaves;j++){
        cumul[j]+=waves->wave[j][i];
        if(cumul[j]>thresh[j]){
          start=i;
          break;
        }
      }
      if(start>=0)break;
    }/*waveform trimming*/
    TIDY(cumul);
    TIDY(thresh);

    start-=buff/gediIO->res;
    if(start<0)start=0;
  }else{
    /*find maximum*/
    ASSIGN_CHECKNULL_RETINT(tot,falloc(1,"tot",0));
    tot[0]=-10000.0;
    for(i=0;i<waves->nBins;i++){
      if(fabs(waves->wave[0][i])>tot[0])tot[0]=fabs(waves->wave[0][i]);
    }
    start=0;
    ASSIGN_CHECKNULL_RETINT(thresh,falloc(1,"thresh",0));
    thresh[0]=tot[0]*0.000001;
    TIDY(tot);
    for(i=0;i<waves->nBins;i++){
      if(fabs(waves->wave[0][i])>thresh[0]){
        start=i-1;
        break;
      }
    }
    if(start<0)start=0;
    buff=0.0;
    TIDY(thresh);
  }

  /*copy data*/
  hdfData->z0[numb]=waves->maxZ-(float)start*gediIO->res;
  hdfData->zN[numb]=hdfData->z0[numb]-(float)hdfData->nBins[0]*gediIO->res;
  if(gediRat->geoCoords){
    hdfData->lon[numb]=gediRat->geoCoords[waveNumb][0];
    hdfData->lat[numb]=gediRat->geoCoords[waveNumb][1];
  }else{
    hdfData->lon[numb]=gediRat->coord[0];
    hdfData->lat[numb]=gediRat->coord[1];
  }
  hdfData->beamDense[numb]=gediRat->beamDense;
  hdfData->pointDense[numb]=gediRat->pointDense;
  hdfData->zen[numb]=waves->meanScanAng;;


  /*ID*/
  if(gediRat->doGrid)snprintf(thisWaveID, sizeof(thisWaveID), "%s.%d.%d", waveID, (int)gediRat->coord[0], (int)gediRat->coord[1]);
  else if(gediRat->waveIDlist)strcpy(thisWaveID,gediRat->waveIDlist[waveNumb]);
  else if(useID)strcpy(thisWaveID,waveID);
  else                  snprintf(thisWaveID, sizeof(thisWaveID), "%d", numb);
  idLength=(hdfData->idLength<((int)strlen(thisWaveID)+1))?hdfData->idLength:(int)strlen(thisWaveID)+1;
  memcpy(&hdfData->waveID[numb*hdfData->idLength],thisWaveID,idLength);

  /*waveform*/
  nBins=(hdfData->nBins[0]<(waves->nBins-start))?hdfData->nBins[0]:waves->nBins-start;
  for(j=0;j<hdfData->nTypeWaves;j++){
    memcpy(&hdfData->wave[j][numb*hdfData->nBins[0]],&waves->wave[j][start],nBins*sizeof(float));
  }

  /*ground variables if using*/
  if(gediIO->ground){
    hdfData->gElev[numb]=waves->gElevSimp;
    hdfData->slope[numb]=waves->gSlopeSimp;
    for(j=0;j<hdfData->nTypeWaves;j++){
      memcpy(&hdfData->ground[j][numb*hdfData->nBins[0]],&waves->ground[j][start],nBins*sizeof(float));
    }
  }

  /*increment counter*/
  (*hdfCount)++;

  return(0);
}/*packGEDIhdf*/


/*##############################################*/
/*set up HDF structure and write header*/

gediHDF *setUpHDF(gediIOstruct *gediIO,gediRatStruct *gediRat,char useID,char *waveID,int *hdfCount,int maxBins)
{
  int i=0;
  gediHDF *hdfData=NULL;

  /*set counter to zero*/
  *hdfCount=0;

  /*allocate space*/
  if(!(hdfData=(gediHDF *)calloc(1,sizeof(gediHDF)))){
    errorf("error control allocation.\n");
    return(NULL);
  }

  /*header*/
  hdfData->nWaves=gediRat->gNx*gediRat->gNy;
  ASSIGN_CHECKNULL_RETNULL(hdfData->nBins, ialloc(1,"nBins",0));
  if(gediIO->pcl==0)hdfData->nBins[0]=(int)((float)maxBins*0.15/gediIO->res);
  else              hdfData->nBins[0]=(int)((gediIO->pulse->x[gediIO->pulse->nBins-1]-gediIO->pulse->x[0])/gediIO->res)*2;
  hdfData->nTypeWaves=gediIO->nTypeWaves;
  hdfData->pSigma=gediIO->pSigma;
  hdfData->fSigma=gediIO->fSigma;

  /*max id label length*/
  if(useID){
    if(gediRat->doGrid){
      hdfData->idLength=(int)strlen(waveID)+1+20;
    }else if(gediRat->waveIDlist){
      hdfData->idLength=-1;
      for(i=0;i<hdfData->nWaves;i++){
        if(((int)strlen(gediRat->waveIDlist[i])+1)>hdfData->idLength)hdfData->idLength=(int)strlen(gediRat->waveIDlist[i])+1;
      }
    }else hdfData->idLength=(int)strlen(waveID)+1;
  }else hdfData->idLength=7;

  /*do we need to record the pulse*/
  //if(gediIO->readPulse){
    hdfData->pRes=gediIO->pRes;
    hdfData->nPbins=gediIO->pulse->nBins;
    ASSIGN_CHECKNULL_RETNULL(hdfData->pulse,falloc((uint64_t)hdfData->nPbins,"hdf pulse",0));
    memcpy(hdfData->pulse,gediIO->pulse->y,sizeof(float)*hdfData->nPbins);
  /*else{
    hdfData->pulse=NULL;
    hdfData->nPbins=0;
  }*/

  /*allocate arrays*/
  ASSIGN_CHECKNULL_RETNULL(hdfData->wave,fFalloc(hdfData->nTypeWaves,"hdf waveforms",0));
  for(i=0;i<hdfData->nTypeWaves;i++) {ASSIGN_CHECKNULL_RETNULL(hdfData->wave[i],falloc((uint64_t)hdfData->nWaves*(uint64_t)hdfData->nBins[0],"hdf waveforms",i+1));}
  if(gediIO->ground){
    ASSIGN_CHECKNULL_RETNULL(hdfData->ground,fFalloc(hdfData->nTypeWaves,"hdf ground waveforms",0));
    for(i=0;i<hdfData->nTypeWaves;i++) {ASSIGN_CHECKNULL_RETNULL(hdfData->ground[i],falloc((uint64_t)hdfData->nWaves*(uint64_t)hdfData->nBins[0],"hdf ground waveforms",i+1));}
  }
  ASSIGN_CHECKNULL_RETNULL(hdfData->z0,falloc((uint64_t)hdfData->nWaves,"hdf z0",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->zN,falloc((uint64_t)hdfData->nWaves,"hdf zN",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->lon,dalloc(hdfData->nWaves,"hdf lon",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->lat,dalloc(hdfData->nWaves,"hdf lat",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->slope,falloc((uint64_t)hdfData->nWaves,"hdf slope",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->gElev,falloc((uint64_t)hdfData->nWaves,"hdf gElev",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->demElev,falloc((uint64_t)hdfData->nWaves,"hdf demElev",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->beamDense,falloc((uint64_t)hdfData->nWaves,"hdf beamDense",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->pointDense,falloc((uint64_t)hdfData->nWaves,"hdf pointDense",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->zen,falloc((uint64_t)hdfData->nWaves,"hdf zen",0));
  ASSIGN_CHECKNULL_RETNULL(hdfData->waveID,challoc(hdfData->nWaves*hdfData->idLength,"hdf waveID",0));

  return(hdfData);
}/*setUpHDF*/


/*##############################################*/
/*update footprint cordinate*/

void updateGediCoord(gediRatStruct *gediRat,int i,int j)
{

  if(gediRat->doGrid){
    gediRat->coord[0]=gediRat->gMinX+(double)i*gediRat->gRes;
    gediRat->coord[1]=gediRat->gMinY+(double)j*gediRat->gRes;
  }else if(gediRat->readALSonce){
    gediRat->coord[0]=gediRat->coords[i][0];
    gediRat->coord[1]=gediRat->coords[i][1];
  }

  return;
}/*updateGediCoord*/


/*####################################*/
/*clear wave structure*/

waveStruct *tidyWaveStruct(waveStruct *waves)
{

  if(waves){
    TTIDY((void **)waves->wave,waves->nWaves);
    TTIDY((void **)waves->canopy,waves->nWaves);
    TTIDY((void **)waves->ground,waves->nWaves);
    TIDY(waves);
  }

  return(waves);
}/*tidyWaveStruct*/


/*####################################*/
/*allocate wave structure*/

waveStruct *allocateGEDIwaves(gediIOstruct *gediIO,gediRatStruct *gediRat,pCloudStruct **data,pointMapStruct *pointmap)
{
  int j=0,numb=0,k=0;
  uint32_t i=0,n=0;
  double maxZ=0,minZ=0;
  double buff=0;
  waveStruct *waves=NULL;
  char hasPoints=0;

  if(!(waves=(waveStruct *)calloc(1,sizeof(waveStruct)))){
    errorf("error waveStruct allocation.\n");
    return(NULL);
  }

  /*determine wave bounds*/
  if(gediIO->pcl==0)buff=35.0;
  else              buff=(double)gediIO->pulse->nBins*(double)gediIO->pRes;
  if(gediIO->pulse||gediIO->pcl)buff+=(double)gediIO->pulse->nBins*(double)gediIO->pRes/2.0;
  minZ=100000000000.0;
  maxZ=-100000000000.0;
  hasPoints=0;

  for(n=0;n<pointmap->nPoints;n++){
    numb=pointmap->fList[n];
    i=pointmap->pList[n];
    if(data[numb]->nPoints>0)hasPoints=1;
    if(data[numb]->z[i]>maxZ)maxZ=data[numb]->z[i];
    if(data[numb]->z[i]<minZ)minZ=data[numb]->z[i];
  }/*bound finding*/

  if(hasPoints==0){
    errorf("No points included\n");
    return(NULL);
  }

  /*determine number of waveform bins*/
  waves->minZ=minZ-buff;
  waves->maxZ=maxZ+buff;
  waves->nBins=(int)((waves->maxZ-waves->minZ)/(double)gediIO->res);

  waves->nWaves=(int)(gediIO->useCount+gediIO->useInt+gediIO->useFrac);
  if(gediRat->readWave)waves->nWaves*=3;  /*if we are using full waveform*/
  ASSIGN_CHECKNULL_RETNULL(waves->wave,fFalloc(waves->nWaves,"result waveform",0));
  for(j=0;j<waves->nWaves;j++){
    ASSIGN_CHECKNULL_RETNULL(waves->wave[j],falloc((uint64_t)waves->nBins,"result waveform",j+1));
    for(k=0;k<waves->nBins;k++)waves->wave[j][k]=0.0;
  }
  if(gediIO->ground){
    ASSIGN_CHECKNULL_RETNULL(waves->canopy,fFalloc(waves->nWaves,"canopy",0));
    ASSIGN_CHECKNULL_RETNULL(waves->ground,fFalloc(waves->nWaves,"ground",0));
    for(j=0;j<waves->nWaves;j++){
      ASSIGN_CHECKNULL_RETNULL(waves->canopy[j],falloc((uint64_t)waves->nBins,"canopy waveform",j+1));
      ASSIGN_CHECKNULL_RETNULL(waves->ground[j],falloc((uint64_t)waves->nBins,"ground waveform",j+1));
      for(k=0;k<waves->nBins;k++)waves->canopy[j][k]=waves->ground[j][k]=0.0;
    }
  }

  return(waves);
}/*allocateGEDIwaves*/


/*################################################################################*/
/*determine ALS coverage*/

void determineALScoverage(gediIOstruct *gediIO,gediRatStruct *gediRat,pCloudStruct **data,pointMapStruct *pointmap)
{
  int i=0;
  int gX=0,gY=0;
  uint32_t j=0,k=0;
  double dx=0,dy=0;
  double sepSq=0;
  float area=0.0;
  char hasPoints=0;


  /*reset counter*/
  for(i=gediRat->gX*gediRat->gX-1;i>=0;i--)gediRat->nGrid[i]=0;

  /*loop over points needed*/
  for(k=0;k<pointmap->nPoints;k++){
    i=pointmap->fList[k];
    j=pointmap->pList[k];
    /*check within bounds*/
    if((data[i]->x[j]>=gediRat->minX)&&(data[i]->x[j]<=gediRat->maxX)&&(data[i]->y[j]>=gediRat->minY)&&(data[i]->y[j]<=gediRat->maxY)){
      dx=data[i]->x[j]-gediRat->coord[0];
      dy=data[i]->y[j]-gediRat->coord[1];
      sepSq=dx*dx+dy*dy;

      /*ground grid for ALS coverage*/
      if(gediRat->normCover||gediRat->checkCover){
        if(data[i]->retNumb[j]==data[i]->nRet[j]){  /*only once per beam*/
          /*mark sampling desnity for normalisation*/
          gX=(int)((data[i]->x[j]-gediRat->g0[0])/(double)gediRat->gridRes);
          gY=(int)((data[i]->y[j]-gediRat->g0[1])/(double)gediRat->gridRes);
          if((gX>=0)&&(gX<gediRat->gX)&&(gY>=0)&&(gY<gediRat->gY)){
            gediRat->nGrid[gY*gediRat->gX+gX]++;
          }
        }
      }/*ground grid for ALS coverage*/

      /*point and beam density*/
      if(sepSq<=gediRat->denseRadSq){
        hasPoints=1;
        gediRat->pointDense+=1.0;
        if(data[i]->retNumb[j]==data[i]->nRet[j])gediRat->beamDense+=1.0;
      }
    }/*bounds check*/
  }/*point loop*/

  area=M_PI*gediRat->denseRadSq;
  gediRat->pointDense/=area;
  gediRat->beamDense/=area;

  if(hasPoints==0){
    for(i=gediRat->gX*gediRat->gX-1;i>=0;i--)gediRat->nGrid[i]=0;
  }

  return;
}/*determineALScoverage*/


/*####################################*/
/*check footprint is covered by ALS*/

void checkFootCovered(gediIOstruct *gediIO,gediRatStruct *gediRat)
{
  int i=0,j=0,nWithin=0;
  int thresh=0,nMissed=0;
  double dX=0,dY=0,sepSq=0;
  double useRad=0,radSq=0;

  useRad=10.0;
  if(useRad>(double)gediIO->fSigma)useRad=(double)gediIO->fSigma;
  radSq=useRad*useRad;

  for(i=0;i<gediRat->gX;i++){
    dX=(double)i*(double)gediRat->gridRes-gediRat->coord[0];
    for(j=0;j<gediRat->gY;j++){
      dY=(double)j*(double)gediRat->gridRes-gediRat->coord[1];
      sepSq=dX*dX+dY*dY;

      if(sepSq<radSq){
        if(gediRat->nGrid[j*gediRat->gX+i]==0)nMissed++;
        nWithin++;
      }
    }/*y loop*/
  }/*x loop*/

  thresh=(int)((float)nWithin*2.0/3.0);
  if(nMissed>thresh){
    errorf("Too many missed %d of %d\n",nMissed,nWithin);
    gediRat->useFootprint=0;
  }else gediRat->useFootprint=1;

  return;
}/*checkFootCovered*/


/*####################################*/
/*set deconvolution parameters for GEDI*/

denPar *setDeconForGEDI(gediRatStruct *gediRat)
{
  void setDenoiseDefault(denPar *);
  int readPulse(denPar *);
  denPar *decon=NULL;

  /*set defaults*/
  if(!(decon=(denPar *)calloc(1,sizeof(denPar)))){
    errorf("error decon structure allocation.\n");
    return(NULL);
  }
  setDenoiseDefault(decon);

  /*particular values for here*/
  decon->deChang=pow(10.0,-8.0);  /*change between decon iterations to stop*/
  decon->thresh=17.0;
  decon->meanN=13.0;
  strcpy(decon->pNamen,"/Users/stevenhancock/data/bess/leica_shape/leicaPulse.dat");  /*pulse filename*/
  decon->deconMeth=0;     /*Gold's method*/
  decon->pScale=1.0;
  decon->noiseTrack=0;

  /*read system pulse*/
  ISINTRETNULL(readPulse(decon));

  return(decon);
}/*setDeconForGEDI*/


/*####################################*/
/*GEDI wave from ALS waveforms*/

int gediFromWaveform(pCloudStruct *data,uint32_t i,float rScale,waveStruct *waves,gediRatStruct *gediRat,gediIOstruct *gediIO)
{
  int j=0,bin=0;
  int buffBins=0;
  uint32_t waveLen=0;
  float grad[3],*smoothed=NULL,*floWave=NULL;
  float *processed=NULL,*smooPro=NULL;
  float *smooth(float,int,float *,float);
  double x=0,y=0,z=0;
  unsigned char *wave=NULL,*temp=NULL;

  for(j=0;j<3;j++)grad[j]=data->grad[j][i];
  ASSIGN_CHECKNULL_RETINT(wave,readLasWave(data->waveMap[i],data->waveLen[i],data->ipoo,data->waveStart));

  /*buffer to give space for smoothing*/
  if(gediIO->pcl==0)buffBins=80;
  else              buffBins=0;
  waveLen=data->waveLen[i]+(uint32_t)(2*buffBins);
  ASSIGN_CHECKNULL_RETINT(temp,uchalloc((uint64_t)waveLen,"temp waveform",0));
  for(j=0;j<buffBins;j++)temp[j]=(unsigned char)gediRat->meanN;
  for(j=0;j<(int)data->waveLen[i];j++)temp[j+buffBins]=wave[j];
  for(j=(int)data->waveLen[i]+buffBins;j<(int)waveLen;j++)temp[j]=(unsigned char)gediRat->meanN;
  TIDY(wave);
  wave=temp;
  temp=NULL;

  /*deconvolve and reconvolve*/
  if(gediRat->indDecon){
    ASSIGN_CHECKNULL_RETINT(processed,processWave(wave,(int)waveLen,gediRat->decon,1.0));
    ASSIGN_CHECKNULL_RETINT(smooPro,smooth(gediIO->pSigma,(int)waveLen,processed,gediIO->res));
  }

  /*convolve with GEDI pulse*/
  ASSIGN_CHECKNULL_RETINT(floWave,falloc((uint64_t)waveLen,"",0));
  for(j=0;j<(int)waveLen;j++)floWave[j]=(float)wave[j]-gediRat->meanN;
  ASSIGN_CHECKNULL_RETINT(smoothed,smooth(gediIO->pSigma,(int)waveLen,floWave,gediIO->res));
  TIDY(floWave);

  /*add up*/
  for(j=0;j<(int)waveLen;j++){
    binPosition(&x,&y,&z,j-buffBins,data->x[i],data->y[i],data->z[i],data->time[i],&(grad[0]));
    bin=(int)((waves->maxZ-z)/(double)gediIO->res);
    if((bin>=0)&&(bin<waves->nBins)){
      waves->wave[3][bin]+=((float)wave[j]-gediRat->meanN)*rScale;  /*with ALS pulse*/
      waves->wave[4][bin]+=smoothed[j]*rScale;
      if(gediRat->doDecon){
        if(gediRat->indDecon){
          waves->wave[5][bin]+=smooPro[j]*rScale;
          waves->wave[6][bin]+=processed[j]*rScale;
        }
        waves->wave[7][bin]+=((float)wave[j]-gediRat->meanN)*rScale;
      }
    }
  }/*wave bin loop*/

  TIDY(wave);
  TIDY(smooPro);
  TIDY(smoothed);
  TIDY(processed);
  return(0);
}/*gediFromWaveform*/


/*################################################################################*/
/*make waveform from point cloud*/

int waveFromPointCloud(gediRatStruct *gediRat, gediIOstruct *gediIO,pCloudStruct **data,waveStruct *waves,pointMapStruct *pointmap)
{
  int numb=0,bin=0,j=0;
  int gX=0,gY=0,n=0;
  int xInd=0,yInd=0;
  uint32_t i=0,k=0;
  double sep=0;
  double dX=0,dY=0;
  double totGround=0;     /*contrbution to ground estimate*/
  float refl=0,rScale=0,fracHit=0,totAng=0;

  /*reset mean scan angle*/
  waves->meanScanAng=totAng=0.0;
  /*ground elevation estimate*/
  if(gediIO->ground){
    waves->gElevSimp=0.0;
    totGround=0.0;
  }

  /*make waves*/
  for(n=0;n<gediRat->nLobes;n++){
    for(k=0;k<pointmap->nPoints;k++){
      numb=pointmap->fList[k];
      i=pointmap->pList[k];

      /*determine laser intensity at this point*/
      dX=data[numb]->x[i]-gediRat->lobe[n].coord[0];
      dY=data[numb]->y[i]-gediRat->lobe[n].coord[1];
      if(gediRat->defWfront==0){    /*symmetric wavefront*/
        sep=sqrt(dX*dX+dY*dY);
        if(gediRat->topHat==0)rScale=(float)gaussian(sep,(double)gediRat->lobe[n].fSigma,0.0);
        else{
          if(sep<=gediRat->lobe[n].maxSepSq)rScale=1.0;
          else                             rScale=0.0;
        }
      }else{     /*read assymmetric footprint*/
        xInd=(int)((dX*cos(gediRat->lobeAng)+dY*sin(gediRat->lobeAng))/(double)gediRat->wavefront->res)+gediRat->wavefront->x0;
        yInd=(int)((dY*cos(gediRat->lobeAng)-dX*sin(gediRat->lobeAng))/(double)gediRat->wavefront->res)+gediRat->wavefront->y0;
        if((xInd>=0)&&(xInd<gediRat->wavefront->nX)&&(yInd>=0)&&(yInd<gediRat->wavefront->nY)){
          if(gediRat->wavefront->front[xInd][yInd]>0.0)rScale=gediRat->wavefront->front[xInd][yInd];
          else rScale=0.0;
        }else rScale=0.0;
      }/*determine laser intensity at this point*/

      if(rScale>gediRat->iThresh){  /*if bright enough to matter*/
        /*scale by sampling density*/
        if(gediRat->normCover){
          gX=(int)((data[numb]->x[i]-gediRat->g0[0])/(double)gediRat->gridRes);
          gY=(int)((data[numb]->y[i]-gediRat->g0[1])/(double)gediRat->gridRes);
          if((gX>=0)&&(gX<gediRat->gX)&&(gY>=0)&&(gY<gediRat->gY)){
            if(gediRat->nGrid[gY*gediRat->gX+gX]>0)rScale/=(float)gediRat->nGrid[gY*gediRat->gX+gX];
          }
        }/*scale by sampling density*/


        /*discrete return*/
        refl=(float)data[numb]->refl[i]*rScale;
        if(data[numb]->nRet[i]>0)fracHit=1.0/(float)data[numb]->nRet[i];
        else                     fracHit=1.0;

        /*if convolving before, smooth now*/
        if(gediRat->pulseAfter==0){
          /*loop over pulse array*/
          for(j=0;j<gediIO->pulse->nBins;j++){
            bin=(int)floor((waves->maxZ-data[numb]->z[i]+(double)gediIO->pulse->x[j])/(double)gediIO->res);
            if((bin>=0)&&(bin<waves->nBins)){
              if(gediIO->useInt)waves->wave[0][bin]+=refl*gediIO->pulse->y[j];
              if(gediIO->useCount)waves->wave[(int)gediIO->useInt][bin]+=rScale*gediIO->pulse->y[j];
              if(gediIO->useFrac)waves->wave[(int)(gediIO->useCount+gediIO->useInt)][bin]+=rScale*fracHit*gediIO->pulse->y[j];
              if(gediIO->ground){
                if(data[numb]->class[i]==2){
                  if(gediIO->useInt)waves->ground[0][bin]+=refl*gediIO->pulse->y[j];
                  if(gediIO->useCount)waves->ground[(int)gediIO->useInt][bin]+=rScale*gediIO->pulse->y[j];
                  if(gediIO->useFrac)waves->ground[(int)(gediIO->useCount+gediIO->useInt)][bin]+=rScale*fracHit*gediIO->pulse->y[j];
                }else{
                  if(gediIO->useInt)waves->canopy[0][bin]+=refl*gediIO->pulse->y[j];
                  if(gediIO->useCount)waves->canopy[(int)gediIO->useInt][bin]+=rScale*gediIO->pulse->y[j];
                  if(gediIO->useFrac)waves->canopy[(int)(gediIO->useCount+gediIO->useInt)][bin]+=rScale*fracHit*gediIO->pulse->y[j];
                }
              }/*ground recording if needed*/
            }/*bin bound check*/
          }/*pulse bin loop*/
        }else{   /*bin up to smooth later*/
          bin=(int)floor((waves->maxZ-data[numb]->z[i])/(double)gediIO->res);
          if((bin>=0)&&(bin<waves->nBins)){
            if(gediIO->useInt)waves->wave[0][bin]+=refl;
            if(gediIO->useCount)waves->wave[(int)gediIO->useInt][bin]+=rScale;
            if(gediIO->useFrac)waves->wave[(int)(gediIO->useCount+gediIO->useInt)][bin]+=rScale*fracHit;
            if(gediIO->ground){
              if(data[numb]->class[i]==2){
                if(gediIO->useInt)waves->ground[0][bin]+=refl;
                if(gediIO->useCount)waves->ground[(int)gediIO->useInt][bin]+=rScale;
                if(gediIO->useFrac)waves->ground[(int)(gediIO->useCount+gediIO->useInt)][bin]+=rScale*fracHit;
              }else{
                if(gediIO->useInt)waves->canopy[0][bin]+=refl;
                if(gediIO->useCount)waves->canopy[(int)gediIO->useInt][bin]+=rScale;
                if(gediIO->useFrac)waves->canopy[(int)(gediIO->useCount+gediIO->useInt)][bin]+=rScale*fracHit;
              }
            }
          }/*bin check*/
        }/*apply pulse before or after*/
        if(gediIO->ground){
          if(data[numb]->class[i]==2){
            waves->gElevSimp+=rScale*data[numb]->z[i];
            totGround+=rScale;
          }
        }
        waves->meanScanAng+=rScale*(float)abs((int)data[numb]->scanAng[i]);
        totAng+=rScale;

        /*full-waveform*/
        if(gediRat->readWave&&data[numb]->hasWave){
          if(data[numb]->packetDes[i]){  /*test for waveform*/
            ISINTRETINT(gediFromWaveform(data[numb],i,rScale,waves,gediRat,gediIO));
          }
        }/*waveform test*/
      }
    }/*point loop*/

    /*if applying pulse after, smooth*/
    if(gediRat->pulseAfter) {ISINTRETINT(applyPulseShape(gediIO,gediRat,waves));}
  }/*lobe loop*/

  /*normalise mean scan angle*/
  if(totAng>0.0)waves->meanScanAng/=totAng;
  if(totGround>=0.0)waves->gElevSimp/=totGround;
  else              waves->gElevSimp=-9999.0;

  return(0);
}/*waveFromPointCloud*/


/*################################################################################*/
/*apply pulse on binned waveform*/

int applyPulseShape(gediIOstruct *gediIO,gediRatStruct *gediRat,waveStruct *waves)
{
  int i=0,j=0,k=0;
  int bin=0;
  int pclSbin=0,pclEbin=0;  /*start and end bins for pcl*/
  float **temp=NULL;
  float **tempGr=NULL;
  float **tempC=NULL;
  float contN=0;
  float minPulse=0;

  /*allocate temporary space*/
  ASSIGN_CHECKNULL_RETINT(temp,fFalloc(waves->nWaves,"temp waves",0));
  if(gediIO->ground){
    ASSIGN_CHECKNULL_RETINT(tempGr,fFalloc(waves->nWaves,"temp ground waves",0));
    ASSIGN_CHECKNULL_RETINT(tempC,fFalloc(waves->nWaves,"temp canopy waves",0));
  }

  /*find min pulse*/
  minPulse=10000.0;
  for(j=0;j<gediIO->pulse->nBins;j++){
    if(gediIO->pulse->y[j]<minPulse)minPulse=gediIO->pulse->y[j];
  }


  /*allocate temporary space*/
  for(k=0;k<waves->nWaves;k++){
    ASSIGN_CHECKNULL_RETINT(temp[k],falloc((uint64_t)waves->nBins,"temp waves",k+1));
    if(gediIO->ground){
      ASSIGN_CHECKNULL_RETINT(tempGr[k],falloc((uint64_t)waves->nBins,"temp ground waves",k+1));
      ASSIGN_CHECKNULL_RETINT(tempC[k],falloc((uint64_t)waves->nBins,"temp canopy waves",k+1));
    }
  }/*allocate temporary space*/

  /*start and end bounds if needed*/
  if(gediIO->pcl){
    pclSbin=(int)floor(((float)gediIO->pulse->nBins*gediIO->pRes)/gediIO->res);
    pclEbin=waves->nBins-(int)floor(((float)gediIO->pulse->nBins*gediIO->pRes)/gediIO->res);
  }

  /*smooth by pulse shape*/
  /*loop over methods*/
  for(k=0;k<waves->nWaves;k++){

    /*loop over waveform bins*/
    for(i=0;i<waves->nBins;i++){
      contN=0.0;    /*reset counters*/
      temp[k][i]=0.0;
      if(gediIO->ground){
        tempGr[k][i]=0.0;
        tempC[k][i]=0.0;
      }

      /*is PCL, only convolve areas that completely overlap*/
      if(gediIO->pcl){
        if((i<pclSbin)||(i>=pclEbin)){
          continue;
        }
      }

      /*loop over pulse bins*/
      for(j=0;j<gediIO->pulse->nBins;j++){

        /*waveform array bin*/
        if(!gediIO->pcl)bin=i+(int)floor((float)(gediIO->pulse->centBin-j)*gediIO->pRes/gediIO->res);
        else            bin=i+(int)floor((float)(j-gediIO->pulse->centBin)*gediIO->pRes/gediIO->res);

        /*are we within the pulse array?*/
        if((bin>=0)&&(bin<waves->nBins)){
          /*add up contributions*/
          temp[k][i]+=waves->wave[k][bin]*gediIO->pulse->y[j];
          if(gediIO->ground){
            tempGr[k][i]+=waves->ground[k][bin]*gediIO->pulse->y[j];
            tempC[k][i]+=waves->canopy[k][bin]*gediIO->pulse->y[j];
          }
          contN+=1.0; //gediIO->pulse->y[j]-minPulse;
        }/*bin bound check*/
      }/*pulse bin loop*/

      /*normalise*/
      if(contN>0.0){
        temp[k][i]/=contN;
        if(gediIO->ground){
          tempGr[k][i]/=contN;
          tempC[k][i]/=contN;
        }
      }/*normalisation step*/
    }/*bin loop*/
  }/*type loop*/

  /*transfer arrays*/
  TTIDY((void **)waves->wave,waves->nWaves);
  TTIDY((void **)waves->ground,waves->nWaves);
  TTIDY((void **)waves->canopy,waves->nWaves);
  waves->wave=temp;
  waves->ground=tempGr;
  waves->canopy=tempC;
  temp=NULL;
  tempC=NULL;
  tempGr=NULL;
  return(0);
}/*applyPulseShape*/


/*################################################################################*/
/*make a map of voxel gaps*/

int voxelGap(gediRatStruct *gediRat,gediIOstruct *gediIO,pCloudStruct **data,waveStruct *waves)
{
  int i=0,vInd=0;
  int xBin=0,yBin=0,zBin=0;
  uint32_t j=0;
  double bounds[6];
  voxStruct *vox=NULL;
  voxStruct *tidyVox(voxStruct *);
  voxStruct *voxAllocate(int,float *,double *,char);
  int countVoxGap(double,double,double,float *,voxStruct *,int,int,float,int);


  bounds[0]=gediRat->minX;
  bounds[1]=gediRat->minY;
  bounds[2]=waves->minZ;
  bounds[3]=gediRat->maxX;
  bounds[4]=gediRat->maxY;
  bounds[5]=waves->maxZ;


  /*first make a voxel map*/
  ASSIGN_CHECKNULL_RETINT(vox,voxAllocate(1,&(gediRat->vRes[0]),&(bounds[0]),0));

  for(i=0;i<gediIO->nFiles;i++){ /*file loop*/
    for(j=0;j<data[i]->nPoints;j++){  /*point loop*/
      ISINTRETINT(countVoxGap(data[i]->x[j],data[i]->y[j],data[i]->z[j],&(data[i]->grad[j][0]),vox,1,1,gediRat->beamRad,0));
    }/*point loop*/
  }/*file loop*/

  /*calculate gap fraction for each return*/
  for(i=0;i<gediIO->nFiles;i++){ /*file loop*/
    ASSIGN_CHECKNULL_RETINT(data[i]->gap,falloc((uint64_t)data[i]->nPoints,"point gaps",i+1));
    for(j=0;j<data[i]->nPoints;j++){  /*point loop*/
      xBin=(int)((data[i]->x[j]-vox->bounds[0])/vox->res[0]+0.5);
      yBin=(int)((data[i]->y[j]-vox->bounds[1])/vox->res[1]+0.5);
      zBin=(int)((data[i]->z[j]-vox->bounds[2])/vox->res[2]+0.5);
      vInd=xBin+vox->nX*yBin+vox->nX*vox->nY*zBin;


      if((vox->hits[0][vInd]+vox->miss[0][vInd])>0.0)data[i]->gap[j]=vox->hits[0][vInd]/(vox->hits[0][vInd]+vox->miss[0][vInd]);
      else                                           data[i]->gap[j]=1.0;
    }/*point loop*/
  }/*file loop*/

  vox=tidyVox(vox);
  return(0);
}/*voxelGap*/


/*################################################################################*/
/*make waveforms accounting for shadowing*/

int waveFromShadows(gediRatStruct *gediRat,gediIOstruct *gediIO,pCloudStruct **data,waveStruct *waves,pointMapStruct *pointmap)
{
  int i=0;
  float **tempWave=NULL;
  //float iRes=0,grad[3];
  float grad[3];
  int voxelGap(gediRatStruct *,gediIOstruct *,pCloudStruct **,waveStruct *);
  rImageStruct *rImage=NULL;    /*range image, a stack nBins long*/
  lidVoxPar lidPar;

  errorf("Silouhette images do not currently work with octrees\n");
  return(-1);

  /*iRes=0.02;*/
  grad[0]=grad[1]=0.0;
  grad[2]=-1.0;

  /*set lidar parameters for a downwards looking ALS*/
  lidPar.minRefl=1.0;             /*minimum refletance value to scale between 0 and 1*/
  lidPar.maxRefl=1.0;             /*maximum refletance value to scale between 0 and 1*/
  lidPar.appRefl=1.0 ;            /*scale between TLS reflectance and size*/
  lidPar.beamTanDiv=0.0;          /*tan of beam divergence*/
  lidPar.beamRad=gediRat->beamRad; /*start radius*/
  lidPar.minGap=0.00001;          /*minimum gap fraction correction to apply*/

  /*gap fraction from voxelising data*/
  ISINTRETINT(voxelGap(gediRat,gediIO,data,waves));

  /*create images*/
  /*rImage=allocateRangeImage(gediIO->nFiles,data,gediIO->pRes*4.0,iRes,&(grad[0]),gediRat->coord[0],gediRat->coord[1],waves->maxZ);*/
  /*rImage=allocateRangeImage(gediIO->nFiles,data,NULL,0.15,0.01,&(grad[0]),gediRat->coord[0],gediRat->coord[1],waves->maxZ,NULL);*/
  errorf("This method is no longer operational. Do not use\n");
  return(-1);


  ISINTRETINT(silhouetteImage(gediIO->nFiles,data,NULL,rImage,&lidPar,NULL,0,NULL));


  /*convert images to waveform*/
  ASSIGN_CHECKNULL_RETINT(tempWave,fFalloc(2,"",0));
  for(i=0;i<2;i++) {ASSIGN_CHECKNULL_RETINT(tempWave[i],falloc((uint64_t)rImage->nBins,"",i+1));}
  waveFromImage(rImage,tempWave,1,gediIO->fSigma);
  for(i=0;i<rImage->nBins;i++)msgf("%f %f %f\n",waves->maxZ-(double)i*rImage->rRes,tempWave[0][i],tempWave[1][i]);
  TTIDY((void **)tempWave,2);
  tempWave=NULL;


  if(rImage){
    TTIDY((void **)rImage->image,rImage->nBins);
    rImage->image=NULL;
    TIDY(rImage);
  }
  return(0);
}/*waveFromShadows*/


/*####################################*/
/*clean outlier points from waveform*/

int cleanOutliers(waveStruct *waves,gediIOstruct *gediIO)
{
  int i=0,j=0,gStart=0;
  char pastGround=0;
  float gGap=0;  /*gap in ground return*/
  float maxGap=0;
  float maxGround=0,gThresh=0;
  float max=0,thresh=0;

  if(!gediIO->ground){
    errorf("No need to clean without ground\n");
    return(-1);
  }

  maxGap=10.0;  /*maximum permittable gap in the ground return*/
  gGap=0.0;
  pastGround=0;

  /*determine max ground and max return*/
  maxGround=max=0.0;
  for(i=0;i<waves->nBins;i++){
    if(waves->ground[(int)gediIO->useInt][i]>maxGround)maxGround=waves->ground[(int)gediIO->useInt][i];  /*note that this uses the count method by default*/
    if(waves->wave[(int)gediIO->useInt][i]>max)max=waves->wave[(int)gediIO->useInt][i];
  }
  gThresh=maxGround*0.01;
  thresh=max*0.001;

  for(i=0;i<waves->nBins;i++){
    if(waves->ground[(int)gediIO->useInt][i]>=gThresh){
      if(pastGround==0)gStart=i;
      pastGround=1;
      gGap=0.0;
    }else{
      if(pastGround)gGap+=gediIO->res;
    }
    if(gGap>maxGap){  /*too big a break, delete*/
      waves->groundBreakElev=waves->maxZ-(double)i*(double)gediIO->res;
      for(;i<waves->nBins;i++){
        waves->canopy[0][i]=waves->canopy[(int)gediIO->useInt][i]=waves->ground[0][i]=waves->ground[(int)gediIO->useInt][i]=0.0;
        for(j=0;j<waves->nWaves;j++)waves->wave[j][i]=0.0;
      }
    }/*too big a break, delete*/
  }

  /*look for above canopy outliers*/
  gGap=0.0;
  maxGap=50.0;
  for(i=gStart;i>=0;i--){
    if(waves->wave[(int)gediIO->useInt][i]>=thresh)gGap=0.0;
    gGap+=gediIO->res;

    if(gGap>=maxGap){
      for(;i>=0;i--){
        waves->canopy[0][i]=waves->canopy[(int)gediIO->useInt][i]=waves->ground[0][i]=waves->ground[(int)gediIO->useInt][i]=0.0;
        for(j=0;j<waves->nWaves;j++)waves->wave[j][i]=0.0;
      }
    }
  }

  return(0);
}/*cleanOutliers*/


/*####################################*/
/*deconvolve aggragated wave*/

int processAggragate(gediRatStruct *gediRat,gediIOstruct *gediIO,waveStruct *waves)
{
  int i=0;
  float *processed=NULL,*smooPro=NULL;
  float *smooth(float,int,float *,float);

  /*Add background noise back*/
  for(i=0;i<waves->nBins;i++)waves->wave[7][i]+=gediRat->meanN;

  /*deconvolve and reconvolve*/
  ASSIGN_CHECKNULL_RETINT(processed,processFloWave(&(waves->wave[7][0]),(int)waves->nBins,gediRat->decon,1.0));
  for(i=0;i<waves->nBins;i++)waves->wave[8][i]=processed[i];
  ASSIGN_CHECKNULL_RETINT(smooPro,smooth(gediIO->pSigma,(int)waves->nBins,processed,gediIO->res));
  TIDY(processed);

  TIDY(waves->wave[7]);
  waves->wave[7]=smooPro;
  smooPro=NULL;

  return(0);
}/*processAggragate*/


/*####################################*/
/*make GEDI waveforms*/

waveStruct *makeGediWaves(gediRatStruct *gediRat,gediIOstruct *gediIO,pCloudStruct **data)
{
  int j=0,k=0;
  float tot=0,minInt=0;
  waveStruct *waves=NULL;
  waveStruct *allocateGEDIwaves(gediIOstruct *,gediRatStruct *,pCloudStruct **,pointMapStruct *);
  int processAggragate(gediRatStruct *,gediIOstruct *,waveStruct *);
  void checkFootCovered(gediIOstruct *,gediRatStruct *);
  int cleanOutliers(waveStruct *,gediIOstruct *);
  int waveFromPointCloud(gediRatStruct *,gediIOstruct *,pCloudStruct **,waveStruct *,pointMapStruct *);
  int waveFromShadows(gediRatStruct *,gediIOstruct *,pCloudStruct **,waveStruct *,pointMapStruct *);
  void determineALScoverage(gediIOstruct *,gediRatStruct *,pCloudStruct **,pointMapStruct *);
  denPar *setDeconForGEDI(gediRatStruct *);
  pointMapStruct *findIntersectingMap(gediRatStruct *,gediIOstruct *,pCloudStruct **);
  pointMapStruct *pointmap=NULL;


  /*determine list of points to use*/
  ASSIGN_CHECKNULL_RETNULL(pointmap,findIntersectingMap(gediRat,gediIO,data));

  /*determine ALS coverage*/
  determineALScoverage(gediIO,gediRat,data,pointmap);

  /*check that whole footprint is covered*/
  if(pointmap->nPoints==0)gediRat->useFootprint=0;
  else if(gediRat->checkCover)checkFootCovered(gediIO,gediRat);
  else                   gediRat->useFootprint=1;

  /*only if it contains data*/
  if(gediRat->useFootprint){
    /*allocate*/
    ASSIGN_CHECKNULL_RETNULL(waves,allocateGEDIwaves(gediIO,gediRat,data,pointmap));

    /*set up denoising if using*/
    if(gediRat->doDecon){
      ASSIGN_CHECKNULL_RETNULL(gediRat->decon,setDeconForGEDI(gediRat));
    }

    /*make waves*/
    if(gediRat->useShadow==0){
      ISINTRETNULL(waveFromPointCloud(gediRat,gediIO,data,waves,pointmap));
    } else                     
      ISINTRETNULL(waveFromShadows(gediRat,gediIO,data,waves,pointmap));

    /*clean outliers if needed*/
    if(gediRat->cleanOut&&(!gediIO->pcl)&&(!gediIO->pclPhoton)){
      ISINTRETNULL(cleanOutliers(waves,gediIO));
    } else waves->groundBreakElev=-100000000.0;

    /*deconvolve aggragated waveform*/
    if(gediRat->doDecon){
      ISINTRETNULL(processAggragate(gediRat,gediIO,waves));
    }

    /*normalise integral*/
    for(k=0;k<waves->nWaves;k++){
      tot=0.0;
      for(j=0;j<waves->nBins;j++)tot+=waves->wave[k][j]*gediIO->res;
      if(tot>0.0){
        for(j=0;j<waves->nBins;j++)waves->wave[k][j]/=tot;
        if(gediIO->ground&&(k<3)){
          for(j=0;j<waves->nBins;j++){
            waves->canopy[k][j]/=tot;
            waves->ground[k][j]/=tot;
          }
        }
      }
    }

    /*tidy arrays*/
    if(gediRat->decon){
      TTIDY((void **)gediRat->decon->pulse,2);
      gediRat->decon->pulse=NULL;
      TIDY(gediRat->decon);
    }

    /*find total integral to check that there is signal*/
    minInt=10000.0;
    for(j=0;j<waves->nBins;j++)if(waves->wave[0][j]<minInt)minInt=waves->wave[0][j];
    tot=0.0;
    for(j=0;j<waves->nBins;j++)tot+=(waves->wave[0][j]-minInt)*gediIO->res;
  }/*contains data*/

  /*check whether empty*/
  if(gediRat->useFootprint){
    if((tot<TOL)||(waves->nBins==0))gediRat->useFootprint=0;
  }
  if(pointmap){
    TIDY(pointmap->fList);
    TIDY(pointmap->pList);
    TIDY(pointmap);
  }

  return(waves);
}/*makeGediWaves*/


/*####################################################*/
/*map points and file indices to use*/

pointMapStruct *findIntersectingMap(gediRatStruct *gediRat,gediIOstruct *gediIO,pCloudStruct **data)
{
  int i=0;
  uint32_t j=0,ind=0;
  pointMapStruct *pointmap=NULL;

  /*search octree or copy all points*/
  if(gediRat->useOctree){
    ASSIGN_CHECKNULL_RETNULL(pointmap,mapFromOctree(gediRat->octList,gediRat->nOct,gediRat->octree,gediRat->minX,gediRat->maxX,gediRat->minY,gediRat->maxY));
  }else{   /*use all points*/
    /*allocate space*/
    if(!(pointmap=(pointMapStruct *)calloc(1,sizeof(pointMapStruct)))){
      errorf("error pointMapStruct allocation.\n");
      return(NULL);
    }
    pointmap->nPoints=0;
    pointmap->fList=NULL;
    pointmap->pList=NULL;

    for(i=0;i<gediIO->nFiles;i++){
      if(data[i]->nPoints==0)continue;
      if(!(pointmap->fList=(int *)realloc(pointmap->fList,(pointmap->nPoints+data[i]->nPoints)*sizeof(int)))){
        errorf("Error allocating memory\n");
        return(NULL);
      }
      if(!(pointmap->pList=(uint32_t *)realloc(pointmap->pList,(pointmap->nPoints+data[i]->nPoints)*sizeof(uint32_t)))){
        errorf("Error allocating memory\n");
        return(NULL);
      }
      for(j=0;j<+data[i]->nPoints;j++){
        ind=pointmap->nPoints+j;
        pointmap->fList[ind]=i;
        pointmap->pList[ind]=j;
      }
      pointmap->nPoints+=data[i]->nPoints;
    }
  }
  return(pointmap);
}/*findIntersectingMap*/


/*####################################################*/
/*allocate space and copy wavefront filename*/

wFrontStruct *copyFrontFilename(char *namen)
{
  wFrontStruct *wavefront=NULL;

  if(!(wavefront=(wFrontStruct *)calloc(1,sizeof(wFrontStruct)))){
    errorf("error in wavefront allocation.\n");
    return(NULL);
  }

  wavefront->front=NULL;
  strcpy(wavefront->frontFile,namen);

  return(wavefront);
}/*copyFrontFilename*/


/*####################################################*/
/*determine true canopy cover*/

float waveformTrueCover(dataStruct *data,gediIOstruct *gediIO,float rhoRatio)
{
  int i=0;
  float totC=0,totG=0;
  float cov=0;
  float minG=0,minC=0;

  if(gediIO->ground){
    totG=totC=0.0;
    minG=minC=100000.0;
    for(i=0;i<data->nBins;i++){
     totG+=data->ground[data->useType][i];
     totC+=data->wave[data->useType][i]-data->ground[data->useType][i];
     if(data->ground[data->useType][i]<minG)minG=data->ground[data->useType][i];
     if(data->wave[data->useType][i]<minC)minC=data->wave[data->useType][i];
    }
    totC-=minC*(float)data->nBins;
    totG-=minG*(float)data->nBins;
    if((totG+totC)>0.0)cov=totC/(totC+totG*rhoRatio);
    else               cov=-1.0;
  }else{
    cov=-1.0;
  }
  return(cov);
}/*waveformTrueCover*/


/*####################################################*/
/*determine Blair sensitivity metric*/

float findBlairSense(dataStruct *data,gediIOstruct *gediIO)
{
  int i=0;
  float gAmp=0,totE=0;
  float sigEff=0,gArea=0;
  float slope=0,tanSlope=0;
  float blairSense=0;
  float meanN=0,stdev=0;
  float notNeeded=0;
  double nNsig=0,nSsig=0;
  float probNoise=0,probMiss=0;
  float *wave=NULL;
  int meanNoiseStats(float *,uint32_t,float *,float *,float *,float,float,int);
  int gaussThresholds(double,double,double,double,double *,double *,noiseThreshStruct *);

  /*set sigma*/
  if(data->fSigma>0.0)gediIO->linkFsig=data->fSigma;

  /*noised if available. Raw if not*/
  if(data->noised)wave=data->noised;
  else            wave=data->wave[data->useType];

  /*determine noise stats for sensitivity metric. Note this is using the threshold to get the mean and stdev*/
  ISINTRETFLT(meanNoiseStats(wave, (uint32_t)data->nBins, &meanN, &stdev, &notNeeded, -1.0, 1.0, (int)(gediIO->den->statsLen / gediIO->res)));
  stdev-=meanN;

  /*total energy*/
  totE=0.0;
  for(i=0;i<data->nBins;i++)totE+=wave[i]-meanN;
  totE*=gediIO->res;
  wave=NULL;

  if(stdev>0.0){
    probNoise=0.05;
    probMiss=0.1;
    ISINTRETFLT(gaussThresholds(1.0, XRES, (double)probNoise, (double)probMiss, &nNsig, &nSsig, &gediIO->noiseSigs));

    slope=2.0*M_PI/180.0;
    tanSlope=sin(slope)/cos(slope);
    gAmp=(float)(nNsig+nSsig)*stdev;
    sigEff=sqrt(gediIO->linkPsig*gediIO->linkPsig+gediIO->linkFsig*gediIO->linkFsig*tanSlope*tanSlope);
    gArea=gAmp*sigEff*sqrt(2.0*M_PI)/totE;

    if(gArea>0.0)(blairSense)=1.0-gArea;
    else         (blairSense)=0.0;
  }else blairSense=1.0;

  return(blairSense);
}/*findBlairSense*/


/*####################################################*/
/*reshape waveform for rhoV rhoG*/

int modifyGroundRho(dataStruct *data,float scaleRhoVrhoG)
{
  int i=0,j=0;
  float newWave=0;

  /*check there is a ground*/
  if(data->ground==NULL){
    errorf("Cannot modify ground rates without ground classification\n");
    return -1;
  }

  /*loop over wave types*/
  for(j=0;j<data->nWaveTypes;j++){
    /*loop over bins*/
    for(i=0;i<data->nBins;i++){
      newWave=(data->wave[j][i]-data->ground[j][i])*scaleRhoVrhoG+data->ground[j][i];
      data->totE[j]+=newWave-data->wave[j][i];
      data->wave[j][i]=newWave;
    }/*bin loop*/
  }/*wave type loop*/

  return 0;
}/*modifyGroundRho*/

/*the end*/
/*####################################################*/

