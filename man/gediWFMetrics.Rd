% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/C_gediWFMetrics.R
\name{gediWFMetrics}
\alias{gediWFMetrics}
\title{GEDI full waveform data processing}
\usage{
gediWFMetrics(
  input,
  outRoot,
  writeFit = FALSE,
  writeGauss = FALSE,
  bounds = NULL,
  ground = FALSE,
  useInt = FALSE,
  useFrac = FALSE,
  rhRes = 5,
  laiRes = 10,
  laiH = 30,
  noRHgauss = FALSE,
  gTol = 0,
  fhdHistRes = 0.001,
  forcePsigma = FALSE,
  bayesGround = FALSE,
  dontTrustGround = FALSE,
  noRoundCoord = FALSE,
  noCanopy = FALSE,
  dcBias = 0,
  nSig = 0,
  hNoise = 0,
  linkNoise = NULL,
  linkFsig = NULL,
  linkPsig = NULL,
  trueSig = NULL,
  bitRate = NULL,
  maxDN = NULL,
  renoise = FALSE,
  newPsig = -1,
  oldPsig = 0.764331,
  addDrift = NULL,
  missGround = FALSE,
  minGap = NULL,
  photonCount = FALSE,
  pcl = FALSE,
  nPhotons = 2.1,
  photonWind = 200,
  noiseMult = 0.1,
  rhoVrhoG = NULL,
  nPhotC = 2.1,
  nPhotG = -1,
  photHDF = FALSE,
  meanN = 0,
  thresh = 1e-14,
  varNoise = FALSE,
  varScale = NULL,
  statsLen = NULL,
  noiseTrack = FALSE,
  sWidth = NULL,
  psWidth = 0,
  msWidth = NULL,
  preMatchF = FALSE,
  postMatchF = FALSE,
  pFile = NULL,
  gWidth = 1.2,
  minGsig = 0.764331,
  minWidth = 0,
  medNoise = FALSE,
  varDrift = NULL,
  driftFac = NULL,
  rhoG = 0.4,
  rhoC = 0.57,
  pSigma = NULL,
  gold = FALSE,
  deconTol = NULL,
  readBinLVIS = FALSE,
  readHDFlvis = FALSE,
  readHDFgedi = TRUE,
  level2 = NULL,
  beamList = NULL,
  skipBeams = NULL,
  readBeams = NULL
)
}
\arguments{
\item{input}{\code{\link[rGEDI:gedi.level1b-class]{rGEDI::gedi.level1b}} (may be a list of objects). Simulated waveform input object(s).}

\item{outRoot}{name. output filename root}

\item{writeFit}{write fitted waveform}

\item{writeGauss}{write Gaussian parameters}

\item{bounds}{minX minY maxX maxY. only analyse data within bounds}

\item{ground}{read true ground from file}

\item{useInt}{use discrete intensity instead of count}

\item{useFrac}{use fractional hits rather than counts}

\item{rhRes}{r. percentage energy resolution of RH metrics}

\item{laiRes}{res. lai profile resolution in metres}

\item{laiH}{h. height to calculate LAI to}

\item{noRHgauss}{do not fit Gaussians}

\item{gTol}{tol. ALS ground tolerance. Used to calculate slope.}

\item{fhdHistRes}{res. waveform intensity resolution to use when calculating FHD from histograms}

\item{forcePsigma}{do not read pulse sigma from file}

\item{bayesGround}{use Bayseian ground finding}

\item{dontTrustGround}{don't trust ground in waveforms, if included}

\item{noRoundCoord}{do not round up coords when outputting}

\item{noCanopy}{do not calculate FHD histograms and LAI profiles}

\item{dcBias}{n. mean noise level}

\item{nSig}{sig. noise sigma}

\item{hNoise}{n. hard threshold noise as a fraction of integral}

\item{linkNoise}{linkM cov. apply Gaussian noise based on link margin at a cover}

\item{linkFsig}{sig. footprint width to use when calculating and applying signal noise}

\item{linkPsig}{sig. pulse width to use when calculating and applying signal noise}

\item{trueSig}{sig. true sigma of background noise}

\item{bitRate}{n. digitisation bit rate}

\item{maxDN}{max. maximum DN}

\item{renoise}{remove noise from truth before applying new noise level}

\item{newPsig}{sig. new value for pulse width, when lengthening pulse}

\item{oldPsig}{sig. old value for pulse width if not defined in waveform file, when lengthening pulse}

\item{addDrift}{xi. apply detector background drift}

\item{missGround}{assume ground is missed to assess RH metrics}

\item{minGap}{gap. delete signal beneath min detectable gap fraction}

\item{photonCount}{output point cloud from photon counting}

\item{pcl}{convert to photon counting pulsecompressed}

\item{nPhotons}{n. mean number of photons}

\item{photonWind}{x. window length for photon counting search, metres}

\item{noiseMult}{x. noise multiplier for photoncounting}

\item{rhoVrhoG}{x. ratio of canopy to ground reflectance at this wavelength. Not different from rhoV and rhoG}

\item{nPhotC}{n. mean number of canopy photons (replaces nPhotons and rhoVrhoG)}

\item{nPhotG}{n. mean number of ground photons (replaces nPhotons and rhoVrhoG)}

\item{photHDF}{write photoncounting}

\item{meanN}{n. mean noise level, if using a predefined mean level}

\item{thresh}{n. noise threshold, if using a predefined noise threshold}

\item{varNoise}{use a variable noise threshold}

\item{varScale}{x. variable noise threshold scale (multiple of stdev above mean to set threshold)}

\item{statsLen}{len. length to calculate noise stats over for varNoise}

\item{noiseTrack}{use noise tracking}

\item{sWidth}{sig. smoothing width, after denoising}

\item{psWidth}{sigma. smoothing width, before denoising}

\item{msWidth}{sig. smoothing width, after noise stats, before denoising}

\item{preMatchF}{matched filter before denoising}

\item{postMatchF}{matched filter after denoising}

\item{pFile}{file. read pulse file, for deconvolution and matched filters}

\item{gWidth}{sig. Gaussian parameter selection smoothing width}

\item{minGsig}{sig. minimum Gaussian sigma to fit}

\item{minWidth}{n. minimum feature width in bins}

\item{medNoise}{use median stats rather than mean}

\item{varDrift}{correct detector drift with variable factor}

\item{driftFac}{xi. fix drift with constant drift factor}

\item{rhoG}{rho. ground reflectance}

\item{rhoC}{rho. canopy reflectance}

\item{pSigma}{sig. pulse width to smooth by if using Gaussian pulse}

\item{gold}{deconvolve with Gold's method}

\item{deconTol}{deconvolution tolerance}

\item{readBinLVIS}{input is an LVIS binary file. Default FALSE.}

\item{readHDFlvis}{read LVIS HDF5 input. Default FALSE.}

\item{readHDFgedi}{read GEDI simulator HDF5 input. Default TRUE.}

\item{level2}{name. level2 filename for LVIS ZG. Default NULL.}

\item{beamList}{character. 0/1 for whether or not to use beams 18, default "11111111"}

\item{skipBeams}{character. list of beam numbers to skip. No spaces between (eg "123")}

\item{readBeams}{character. list of beam numbers to read. No spaces between (eg "123")}
}
\value{
Returns a list of metrics derived from the simulated full waveform.
A text file (txt) containing the metrics will be saved in the output folder (outRoot).
Please see the details section for checking the definition of the metrics.
}
\description{
GEDI full waveform data processing and metrics extraction
}
\details{
a) Metrics descriptions

a.1) Metrics available to GEDI
\itemize{
\item \emph{gHeight} Ground elevation (m) from Gaussian fitting
\item \emph{maxGround} Ground elevation (m) from lowest maximum
\item \emph{inflGround} Ground elevation (m) from inflection points.
\item \emph{signal top} Elevation of first point above noise (may include noise tracking).
\item \emph{signal bottom} Elevation of last return above noise (may include noise tracking).
\item \emph{cover} Canopy cover (fraction) from area of Gaussian fitted ground. Uses rho_v=0.57 and rho_g=0.4.
\item \emph{leading edge ext} Leading edge extent (m), from Lefksy et al (2007).
\item \emph{trailing edge extent} Trailing edge extent (m), from Lefksy et al (2007).
\item \emph{rhGauss 0-100} RH metrics, 0\%-100\%, using ground from Gaussian fitting (m).
\item \emph{rhMax 0-100} RH metrics, 0\%-100\%, using ground from lowest maximum (m).
\item \emph{rhInfl 0-100} RH metrics, 0\%-100\%, using ground from inflection points (m).
\item \emph{gaussHalfCov} Canopy cover (fraction) from double the energy beneath the Gaussian ground. Uses rho_v=0.57 and rho_g=0.4.
\item \emph{maxHalfCov} Canopy cover (fraction) from double the energy beneath the lowest maximum ground. Uses rho_v=0.57 and rho_g=0.4.
\item \emph{infHalfCov} Canopy cover (fraction) from double the energy beneath the inflection point ground. Uses rho_v=0.57 and rho_g=0.4.
\item \emph{bayHalfCov} Canopy cover (fraction) from double the energy beneath the experimental "Bayesian" ground. Uses rho_v=0.57 and rho_g=0.4.
\item \emph{lon} Footprint centre longitude in projection of ALS data (m).
\item \emph{lat} Footprint centre latitude in projection of ALS data (m).
\item \emph{waveEnergy} Total energy within waveform (will be 1 scaled by noise for simulations).
\item \emph{blairSense} Blair's sensitivity metric. Canopy cover at which this SNR would have90\% chance of detecting ground (does not account for rho_v/rho_g).
\item \emph{FHD} Foliage height diversity
\item \emph{niM2} Wenge Ni's biomass metric, equal to the sum of the RH metrics to the power of 2 (unpublished)
\item \emph{niM2.1} Wenge Ni's biomass metric, equal to the sum of the RH metrics to the power of 2.1 (unpublished)
}

a.2) Metrics unavailable to GEDI
\itemize{
\item \emph{wave ID} Waveform label, relates to plot name and footprint number.
\item \emph{true ground} Ground elevation (m) from ALS. Centre of gravity of ground points within footprint
\item \emph{true top} Levation of highest point of waveform (m), without noise. Includes pulse blurring.
\item \emph{ground slope} Effective ground slope (degrees), from width of ground return. Includes roughness.
\item \emph{ALS cover} Canopy cover (fraction) from ALS data. Uses rho_v=0.57 and rho_g=0.4.
\item \emph{rhReal 0-100} RH metrics, 0\%-100\%, using "true" ground from ALS data (m).
\item \emph{groundOverlap} Fraction of ground return overlapping with canopy return. A measure of understorey.
\item \emph{groundMin} Depth of minimum between ground and canopy return. A measure of understorey.
\item \emph{groundInfl} d2y/dx2 of inflection point between ground and canopy return. A measure of understorey.
\item \emph{pointDense} Average ALS point density within GEDI footprint.
\item \emph{beamDense} Average ALS beam density within GEDI footprint.
}

a.3) System settings
\itemize{
\item \emph{pSigma} GEDI system pulse width, sigma (m).
\item \emph{fSigma} GEDI footprint width, sigma (m).
\item \emph{linkM} Link margin if noise is added (db).
\item \emph{linkCov} Canopy cover at which the above link margin is true (fraction).
\item \emph{filename} Name of input waveform file.
}

b) Signal processing description
\itemize{
\item \emph{Gaussian fitting} Used for "gHeight", "rhGauss" and "gaussHalfCov".
The waveform is denoised (mean+5\emph{sigma, noise tracking to avoid truncation), smoothed (pSigma}0.75) and Gaussians fitted with Levenberg-Marquardt optimisation.
The center of the lowest Gaussian containing at least 0.5\% of the waveform energy is selected as the ground.

\item \emph{Maximum} Used for "maxGround", "rhMax" and "maxHalfCov".
The waveform is denoised (mean+5\emph{sigma, noise tracking to avoid truncation), smoothed (pSigma}0.75).
The lowest maximum is taken as the ground.

\item \emph{Inflection points} Used for "inflGround", "rhInfl" and "inflHalfCov".
The waveform is denoised (mean+5\emph{sigma, noise tracking to avoid truncation), smoothed (pSigma}0.75).
The centre of gravity between the lowest two inflection points is taken as the ground.

\item \emph{Half covers} Used for "gaussHalfCov", "maxHalfCov" and "inflHalfCov".
Sum energy beneath estimated ground position.
Double that is the ground energy.
Calculate canopy cover, correcting for rho_v and rho_g.

\eqn{cover = \frac{E_{can}}{E_{can} + E_g*\frac{rho_v}{rho_g}}}

Where Ecan is the canopy energy, Eg is the ground energy, rho_v is the vegetation reflectance and rho_g is the ground reflectance.

\item \emph{Edge extents} These are described in: Lefsky, Michael A., Michael Keller, Yong Pang, Plinio B. De Camargo, and Maria O. Hunter.
"Revised method for forest canopy height estimation from Geoscience Laser Altimeter System waveforms." Journal of Applied Remote Sensing 1, no. 1 (2007): 013537-013537.
}
}
\examples{
\dontshow{
rm(list = ls())
}

libsAvailable <- require(lidR) && require(plot3D)
if (libsAvailable) {
  outdir <- tempdir()

  # Specifying the path to ALS data (zip)
  alsfile_Amazon_zip <- system.file("extdata", "Amazon.zip", package = "rGEDIsimulator")
  alsfile_Savanna_zip <- system.file("extdata", "Savanna.zip", package = "rGEDIsimulator")

  # Unzipping ALS data
  alsfile_Amazon_filepath <- unzip(alsfile_Amazon_zip, exdir = outdir)
  alsfile_Savanna_filepath <- unzip(alsfile_Savanna_zip, exdir = outdir)

  # Reading and plot ALS file (las file)
  als_Amazon <- readLAS(alsfile_Amazon_filepath)
  als_Savanna <- readLAS(alsfile_Savanna_filepath)

  # Extracting plot center geolocations
  xcenter_Amazon <- mean(st_bbox(als_Amazon)[1, ])
  ycenter_Amazon <- mean(st_bbox(als_Amazon)[2, ])
  xcenter_Savanna <- mean(st_bbox(als_Savanna)[1, ])
  ycenter_Savanna <- mean(st_bbox(als_Savanna)[2, ])

  # Simulating GEDI full waveform
  wf_Amazon <- gediWFSimulator(
    input = alsfile_Amazon_filepath,
    output = file.path(outdir, "gediWF_amazon_simulation.h5"),
    coords = c(xcenter_Amazon, ycenter_Amazon)
  )

  wf_Savanna <- gediWFSimulator(
    input = alsfile_Savanna_filepath,
    output = file.path(outdir, "gediWF_Savanna_simulation.h5"),
    coords = c(xcenter_Savanna, ycenter_Savanna)
  )

  # Extracting GEDI full waveform derived metrics without adding noise to the full waveform
  wf_amazon_metrics <- gediWFMetrics(input = wf_Amazon, outRoot = file.path(outdir, "amazon"))
  wf_savanna_metrics <- gediWFMetrics(input = wf_Savanna, outRoot = file.path(outdir, "savanna"))

  metrics <- rbind(wf_amazon_metrics, wf_savanna_metrics)
  rownames(metrics) <- c("Amazon", "Savanna")
  head(metrics)

  # Extracting GEDI full waveform derived metrics after adding noise to the waveform
  wf_amazon_metrics_noise <- gediWFMetrics(
    input = wf_Amazon,
    outRoot = file.path(outdir, "amazon"),
    linkNoise = c(3.0103, 0.95),
    maxDN = 4096,
    sWidth = 0.5,
    varScale = 3
  )

  wf_savanna_metrics_noise <- gediWFMetrics(
    input = wf_Savanna,
    outRoot = file.path(outdir, "savanna"),
    linkNoise = c(3.0103, 0.95),
    maxDN = 4096,
    sWidth = 0.5,
    varScale = 3
  )

  close(wf_Amazon)
  close(wf_Savanna)

  metrics_noise <- rbind(wf_amazon_metrics_noise, wf_savanna_metrics_noise)
  rownames(metrics_noise) <- c("Amazon", "Savanna")
  head(metrics_noise)
}
}
\seealso{
i) Hancock, S., Armston, J., Hofton, M., Sun, X., Tang, H., Duncanson, L.I., Kellner,
J.R. and Dubayah, R., 2019. The GEDI simulator: A large-footprint waveform lidar simulator
for calibration and validation of spaceborne missions. Earth and Space Science.
\doi{10.1029/2018EA000506}

ii) gediSimulator: \url{https://bitbucket.org/StevenHancock/gedisimulator/src/master/}
}
