#ifndef _22E2D91D_B134_3B0A_133E_DE2D88E57819
#define _22E2D91D_B134_3B0A_133E_DE2D88E57819
/*##############################*/
/*# Structures for storing     #*/
/*# simulated GEDI waveforms   #*/
/*# 2017 svenhancock@gmail.com #*/
/*##############################*/

/*#######################################*/
/*# Copyright 2015-2017, Steven Hancock #*/
/*# The program is distributed under    #*/
/*# the terms of the GNU General Public #*/
/*# License.    svenhancock@gmail.com   #*/
/*#######################################*/


/*########################################################################*/
/*# This file is part of the NASA GEDI simulator, gediRat.               #*/
/*#                                                                      #*/
/*# gediRat is free software: you can redistribute it and/or modify      #*/
/*# it under the terms of the GNU General Public License as published by #*/
/*# the Free Software Foundation, either version 3 of the License, or    #*/
/*#  (at your option) any later version.                                 #*/
/*#                                                                      #*/
/*# gediRat is distributed in the hope that it will be useful,           #*/
/*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
/*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
/*#   GNU General Public License for more details.                       #*/
/*#                                                                      #*/
/*#    You should have received a copy of the GNU General Public License #*/
/*#    along with gediRat.  If not, see <http://www.gnu.org/licenses/>.  #*/
/*########################################################################*/

#ifdef USEPHOTON
#include "photonCount.h"
#endif
#include "gediNoise.h"
#include <libLasProcess.h>

/*###########################################################*/
/*LVIS level2 data*/

typedef struct{
  uint64_t numb;      /*number of records*/
  uint32_t *lfid;     /*LVIS file identifier*/
  uint32_t *shotN;    /*LVIS shotnumber*/
  float *zG;          /*ground elevation*/
}lvisL2struct;


/*####################################*/
/*empty structure if photon counting not provided*/

#ifndef USEPHOTON
typedef struct{
  void *nothing;
}photonStruct;
#endif


/*####################################*/
/*structure to hold SNR*/

typedef struct{
  int nWaves;
  /*false rates*/
  float falsePosRate;
  float falseNegRate;
  /*smoothing widths*/
  int nSig;
  float minSig;
  float maxSig;
  float dSig;
  /*min widths*/
  int nMinWid;
  int minWid;
  int maxWid;
  int dWid;
  /*array of SNRs*/
  float ***linkM;   /*link margin per sWidth, per minWidth, per wave*/
  float ***bSense;  /*beam sensitivity per sWidth, per minWidth, per wave*/
  float *cov;       /*canopy cover per wave*/
  float *gWidth;    /*ground width per wave*/
}snrStruct;


/*####################################*/
/*control structure*/

typedef struct{
  /*input/output*/
  gediIOstruct gediIO; /*input/output structure*/
  char outRoot[200];
  FILE *opooGauss;  /*Gaussian parameter output*/
  FILE *opooMet;    /*waveform metric output*/
  int maxGauss;     /*maximum number of Gaussians for output*/

  /*level2 LVIS for ZG*/
  char l2namen[200]; /*list of level2 filenames*/
  char readL2;      /*switch to read L2 or not*/

  /*switches*/
  char writeFit;    /*write fitted wave switch*/
  float rhRes;      /*rh resolution*/
  char bayesGround; /*Bayseian ground finding*/
  char noRHgauss;   /*do not do Gaussian fitting*/
  char renoiseWave; /*remove noise before adding*/
  char readBinLVIS;  /*read binary LVIS rather than a list of ASCII files*/
  char readHDFlvis;  /*read HDF5 LVIS rather than ASCII*/
  char readHDFgedi;  /*read HDF5 GEDI rather than ASCII*/
  char coord2dp;     /*round up coords to 2dp when writing*/
  char useBounds;    /*when we will process only a subset of bounds*/
  char writeGauss;   /*write Gaussian parameters*/
  char noCanopy;     /*output the FHD and LAI profile switch*/
  char readPulse;    /*read pulse from an ASCII file*/
  float laiRes;      /*LAI profile resolution*/
  float maxLAIh;     /*maximum height bin of LAI profile. Put all above this in top bin*/
  char rhNoGround;   /*do not use ground in RH metrics*/
  char onlySNR;      /*only calculate the SNR*/

  /*noise parameters*/
  noisePar noise;  /*noise adding structure*/
  float bThresh;   /*bounds threshold*/

  /*LVIS or HDF data*/
  lvisLGWstruct lvis;   /*LVIS lgw structure*/
  lvisHDF *hdfLvis;     /*LVIS HDF5 structure*/
  lvisL2struct *lvisL2; /*LVIS level2 data*/
  gediHDF *hdfGedi;     /*GEDI HDF5 structure*/

  /*bounds for subsets*/
  double minX;
  double maxX;
  double minY;
  double maxY;

  /*photon counting*/
  char ice2;         /*ICESat-2 mode. GEDI by default*/
  photonStruct photonCount;  /*photon counting structure*/

  /*SNR structure*/
  snrStruct *snr;      /*structure to hold SNR*/

  /*others*/
  float rhoRatio;      /*ratio of canopy to ground reflectance, used only for true canopy cover*/
  float scaleRhoVrhoG; /*this is used to rescale the true waveform*/
  char changeGrRho;    /*switch to rescale the true waveform*/
  float gTol;          /*toleranve used to label ALS ground finding*/
  float zen;           /*zenith angle*/
  float fhdHistRes;    /*resolution for FHD histogram method*/
}control;


/*###########################################################*/
/*Bayseian ground structure*/

typedef struct{
  double gHeight;   /*ground elevation*/
  float cov;        /*canopy cover*/
  float slope;      /*slope, degrees*/
}bGround;


/*###########################################################*/
/*metric structure*/

typedef struct{
  float *rh;        /*rh metrics using Gaussian ground*/
  float *rhMax;     /*rh metrics using max ground*/
  float *rhInfl;    /*rh metrics using inflection ground*/
  float *rhReal;    /*rh metric from real ground*/
  int nRH;          /*number of RH metrics*/
  float FHD;        /*foliage height diversity, all waveform, wave*/
  float FHDhist;    /*foliage height diversity, all waveform, hist*/
  float FHDcan;     /*foliage height diversity, canopy, wave*/
  float FHDcanH;    /*foliage height diversity, canopy, hist*/
  float FHDcanGauss;/*foliage height diversity, canopy from Gaussian fitting, wave*/
  float FHDcanGhist;/*foliage height diversity, canopy from Gaussian fitting, hist*/
  int nLm;          /*number of L-moments*/
  //float *LmomGau;   /*L-moments from Gaussian fit*/
  //float *LmomRea;   /*L-moments from ALS ground*/
  //float *LmomInf;   /*L-moments from inflection point*/
  //float *LmomMax;   /*L-moments from maximum*/
  float cov;        /*canopy cover for gaussian fitting*/
  double gHeight;   /*ground height from Gaussians*/
  float gSlope;     /*slope estimate from Gaussian fitting*/
  double maxGround; /*ground height from maximum*/
  double inflGround;/*ground height from inflection*/
  double tElev;     /*top elevation*/
  double bElev;     /*bottom elevation*/
  float leExt;      /*Lefsky's leading edge extent*/
  float teExt;      /*Lefsky's trailing edge extent*/
  float covHalfG;   /*cover from Bryan's half, Gaussian*/
  float covHalfI;   /*cover from Bryan's half, Inflection*/
  float covHalfM;   /*cover from Bryan's half, maximum*/
  float covHalfB;   /*cover from Bryan's half, Bayesian*/
  float totE;       /*total energy after denoising*/
  float blairSense; /*Blair sensitivity metric*/
  float niM2;       /*Ni metric with c=2*/
  float niM21;      /*Ni metric with c=2.1*/
  float *tLAI;      /*true LAI profile*/
  float *gLAI;      /*LAI profile with Gaussian ground removal*/
  float *hgLAI;     /*LAI profile with halp width ground removal, Gaussian elevation*/
  float *hiLAI;     /*LAI profile with halp width ground removal, inflection elevation*/
  float *hmLAI;     /*LAI profile with halp width ground removal, maximum elevation*/
  int laiBins;      /*number of LAI bins*/

  int nBgr;         /*number of ground estimates*/
  bGround *bGr;     /*Bayesian ground structure*/
  double bayGround; /*Bayesian ground elevation*/
}metStruct;

float *findLAIprofile(float *,float,int,float,int *,double,float,double *,float,float);
char checkUsable(float *,int);
int calculateSNR(control*, dataStruct*, int);
int writeSNR(char*, snrStruct*);
void tidySNR(control*);
int allocateSNR(control *);

#endif /* _22E2D91D_B134_3B0A_133E_DE2D88E57819 */
