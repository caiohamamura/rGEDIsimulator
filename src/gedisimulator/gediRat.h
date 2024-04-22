#ifndef _B0ABC6EC_B134_3B09_11AC_315F983A6CD0
#define _B0ABC6EC_B134_3B09_11AC_315F983A6CD0
/*####################################*/
/*control structure*/

typedef struct
{
  char **inList;
  char outNamen[1000];
  char waveNamen[400];

  /*IO structure*/
  gediIOstruct gediIO;   /*generic IO options*/
  gediRatStruct gediRat; /*simulator options*/

  /*options*/
  char listFiles;     /*list waves only*/
  char overWrite;     /*overwrite old waveform switch*/
  uint64_t pBuffSize; /*point buffer rading size in bytes*/
  char waveID[200];   /*wave ID if we are to use it*/
  char useID;         /*use wave ID*/
  char polyGr;        /*fit a polynomial to the ground*/
  char nnGr;          /*ground DEM from nearest neighbour*/

  /*HDF5 output*/
  char writeHDF; /*write output as hdf5*/
  char writeL1B; /*write L1B HDF5 output format*/
  int maxBins;   /*bins per wave for HDF5 output*/
  int hdfCount;  /*count used footprints*/
} control;

pCloudStruct *readAsciiData(char *inNamen);
#endif /* _B0ABC6EC_B134_3B09_11AC_315F983A6CD0 */
