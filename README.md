![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig1.png)<br/>

[![R-CMD-check](https://github.com/caiohamamura/rGEDIsimulator/actions/workflows/build.yaml/badge.svg?branch=main)](https://github.com/caiohamamura/rGEDIsimulator/actions/workflows/build.yaml)
[![CRAN](https://www.r-pkg.org/badges/version/rGEDIsimulator)](https://cran.r-project.org/package=rGEDIsimulator)
![Github](https://img.shields.io/badge/Github-0.3.0-green.svg)
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rGEDIsimulator)

**rGEDISimulator: GEDI simulator extension for rGEDI**

Authors: Caio Hamamura, Carlos Alberto Silva, Ruben Valbuena, Steven Hancock, Adrian Cardil, Eben N. Broadbent, Danilo R. A. de Almeida, Celso H. L. Silva Junior and Carine Klauberg  

The rGEDIsimulator will provide fullwaveform GEDI data simulation and calculates metrics based on aerial lidar systems data (ALS).

# Getting Started

## Installation

First we need to install rGEDI itself:
```r
# Install from github
install.packages('rGEDI', repos = c('https://carlos-alberto-silva.r-universe.dev', 'https://cloud.r-project.org'))

library(rGEDI)
```    

Then we can install the rGEDI simulator

```r
# Install from github
install.packages('rGEDIsimulator', repos = c('https://caiohamamura.r-universe.dev', 'https://cloud.r-project.org'))

library(rGEDIsimulator)
```    

## Simulating GEDI full-waveform data from Airborne Laser Scanning (ALS) 3-D point cloud and extracting canopy derived metrics
```r
outdir=getwd()

#######
# Herein, we are using only a GEDI sample dataset for this tutorial.
#######
# downloading zip file
download.file("https://github.com/carlos-alberto-silva/rGEDI/releases/download/datasets/examples.zip",destfile=file.path(outdir, "examples.zip"))

# unzip file 
unzip(file.path(outdir,"examples.zip"))


# Specifying the path to ALS data
lasfile_amazon <- file.path(outdir, "Amazon.las")
lasfile_savanna <- file.path(outdir, "Savanna.las")

# Reading and plot ALS file
library(lidR)
library(plot3D)
las_amazon<-readLAS(lasfile_amazon)
las_savanna<-readLAS(lasfile_savanna)

# Extracting plot center geolocations
xcenter_amazon = mean(st_bbox(las_amazon)[c(1, 3)])
ycenter_amazon = mean(st_bbox(las_amazon)[c(2, 4)])
xcenter_savanna = mean(st_bbox(las_savanna)[c(1, 3)])
ycenter_savanna = mean(st_bbox(las_savanna)[c(2, 4)])

# Simulating GEDI full-waveform
wf_amazon<-gediWFSimulator(input=lasfile_amazon,output=file.path(getwd(),"gediWF_amazon_simulation.h5"),coords = c(xcenter_amazon, ycenter_amazon))
wf_savanna<-gediWFSimulator(input=lasfile_savanna,output=file.path(getwd(),"gediWF_savanna_simulation.h5"),coords = c(xcenter_savanna, ycenter_savanna))

# Plotting ALS and GEDI simulated full-waveform
png("gediWf.png", width = 8, height = 6, units = 'in', res = 300)

par(mfrow=c(2,2), mar=c(4,4,0,0), oma=c(0,0,1,1),cex.axis = 1.2)
scatter3D(las_amazon@data$X,las_amazon@data$Y,las_amazon@data$Z,pch = 16,colkey = FALSE, main="",
          cex = 0.5,bty = "u",col.panel ="gray90",phi = 30,alpha=1,theta=45,
          col.grid = "gray50", xlab="UTM Easting (m)", ylab="UTM Northing (m)", zlab="Elevation (m)")

# Simulated waveforms shot_number is incremental beggining from 0
shot_number = 0
simulated_waveform_amazon = getLevel1BWF(wf_amazon, shot_number)
plot(simulated_waveform_amazon, relative=TRUE, polygon=TRUE, type="l", lwd=2, col="forestgreen",
     xlab="", ylab="Elevation (m)", ylim=c(90,140))
grid()
scatter3D(las_savanna@data$X,las_savanna@data$Y,las_savanna@data$Z,pch = 16,colkey = FALSE, main="",
          cex = 0.5,bty = "u",col.panel ="gray90",phi = 30,alpha=1,theta=45,
          col.grid = "gray50", xlab="UTM Easting (m)", ylab="UTM Northing (m)", zlab="Elevation (m)")

shot_number = 0
simulated_waveform_savanna = getLevel1BWF(wf_savanna, shot_number)
plot(simulated_waveform_savanna, relative=TRUE, polygon=TRUE, type="l", lwd=2, col="green",
xlab="Waveform Amplitude (%)", ylab="Elevation (m)", ylim=c(815,835))
grid()
dev.off()
```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig7.png)

## Extracting GEDI full-waveform derived metrics without adding noise to the full-waveform
```
wf_amazon_metrics<-gediWFMetrics(input=wf_amazon,
                                outRoot=file.path(getwd(), "amazon"))
wf_savanna_metrics<-gediWFMetrics(input=wf_savanna,
                                outRoot=file.path(getwd(), "savanna"))

metrics<-rbind(wf_amazon_metrics,wf_savanna_metrics)
rownames(metrics)<-c("Amazon","Savanna")
head(metrics[,1:8])

#                wave ID true ground true top ground slope ALS cover gHeight maxGround inflGround
#Amazon  gedi.BEAM0000.0      -1e+06   133.25       -1e+06        -1   94.93     99.95      95.16
#Savanna gedi.BEAM0000.0      -1e+06   831.47       -1e+06        -1  822.18    822.17     822.25
```
## Extracting GEDI full-waveform derived metrics after adding noise to the full-waveform
```
wf_amazon_metrics_noise<-gediWFMetrics(input=wf_amazon,
                         outRoot=file.path(getwd(), "amazon"),
                         linkNoise= c(3.0103,0.95),
                         maxDN= 4096,
                         sWidth= 0.5,
                         varScale= 3)

wf_savanna_metrics_noise<-gediWFMetrics(
                        input=wf_savanna,
                        outRoot=file.path(getwd(), "savanna"),
                        linkNoise= c(3.0103,0.95),
                        maxDN= 4096,
                        sWidth= 0.5,
                        varScale= 3)

metrics_noise<-rbind(wf_amazon_metrics_noise,wf_savanna_metrics_noise)
rownames(metrics_noise)<-c("Amazon","Savanna")
head(metrics_noise[,1:8])

#         #wave ID true ground true top ground slope ALS cover gHeight maxGround inflGround
# Amazon         0      -1e+06   133.29       -1e+06        -1   99.17     99.99      95.39
# Savanna        0      -1e+06   831.36       -1e+06        -1  822.15    822.21     822.18

```

## Always close gedi objects, so HDF5 files won't be blocked!
```{r cleanup, echo=TRUE, results="hide", error=TRUE}
close(wf_amazon)
close(wf_savanna)
close(gedilevel1b)
close(gedilevel2a)
close(gedilevel2b)
```


# References
Dubayah, R., Blair, J.B., Goetz, S., Fatoyinbo, L., Hansen, M., Healey, S., Hofton, M., Hurtt, G.,         Kellner, J., Luthcke, S., & Armston, J. (2020) The Global Ecosystem Dynamics Investigation:         High-resolution laser ranging of the Earth’s forests and topography. Science of Remote             Sensing, p.100002. https://doi.org/10.1016/j.srs.2020.100002

Hancock, S., Armston, J., Hofton, M., Sun, X., Tang, H., Duncanson, L.I., Kellner,
       J.R. and Dubayah, R., 2019. The GEDI simulator: A large-footprint waveform lidar simulator
       for calibration and validation of spaceborne missions. Earth and Space Science.
       https://doi.org/10.1029/2018EA000506

Silva, C. A.; Saatchi, S.; Alonso, M. G. ; Labriere, N. ; Klauberg, C. ; Ferraz, A. ; Meyer, V. ;        Jeffery, K. J. ; Abernethy, K. ; White, L. ; Zhao, K. ; Lewis, S. L. ; Hudak, A. T. (2018)         Comparison of Small- and Large-Footprint Lidar Characterization of Tropical Forest                 Aboveground Structure and Biomass: A Case Study from Central Gabon. IEEE Journal of Selected       Topics in Applied Earth Observations and Remote Sensing, p. 1-15. https://doi.org/10.1109/JSTARS.2018.2816962

GEDI webpage. Accessed on September 23 2025 https://gedi.umd.edu/   
GEDI L1B Geolocated Waveform Data Global Footprint Level V002. Accessed on September 23 2025 https://doi.org/10.5067/GEDI/GEDI01_B.002
GEDI L2A Elevation and Height Metrics Data Global Footprint Level V002. Accessed on September 23 2025 https://doi.org/10.5067/GEDI/GEDI02_A.002
GEDI L2B Canopy Cover and Vertical Profile Metrics Data Global Footprint Level V002. Accessed on September 23 2025 https://doi.org/10.5067/GEDI/GEDI02_B.002

# Acknowledgements
The University of Maryland and NASA's Goddard Space Flight Center for developing GEDI mission.

We gratefully acknowledge funding from NASA’s Carbon Monitoring Systems, grant NNH15ZDA001N-CMS. Project entitled "Future Mission Fusion for High Biomass Forest Carbon Accounting" led by Dr. Laura Duncanson (lduncans@umd.edu, University of Maryland) and Dr. Lola Fatoyinbo (lola.fatoyinbo@nasa.gov, NASA's Goddard Space Flight Center).

The Brazilian National Council for Scientific and Technological Development (CNPq) for funding the project entitled "Mapping fuel load and simulation of fire behaviour and spread in the Cerrado biome using modeling and remote sensing technologies" and leaded by Prof. Dr. Carine Klauberg (carine_klauberg@hotmail.com) and Dr. Carlos Alberto Silva
(carlos_engflorestal@outlook.com).

# Getting Help
The best place to get help from community is StackExchange:
<https://gis.stackexchange.com/questions/tagged/gedi>. 
Before posting there, make sure your question hasn't already been answered.
Also don't forget to add relevant tags such as `gedi` and `rgedi`.

# Citing rGEDIsimulator
Hamamura,C.; Silva,C.A; Hancock,S.; Valbuena, R.; Cardil,A.; Broadbent, E. N.; Almeida,D.R.A.; Silva Junior, C.H.L; Klauberg, C. NASA's Global Ecosystem Dynamics Investigation (GEDI) Simulator for ALS Data.
version 0.3.0, accessed on April. 04 2024, available at: <https://CRAN.R-project.org/package=rGEDIsimulator>

# Disclaimer
**rGEDI package has not been developted by the GEDI team. It comes with no guarantee, expressed or implied, and the authors hold no responsibility for its use or reliability of its outputs.**

