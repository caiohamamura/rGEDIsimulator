#' GEDI full waveform data simulation
#'
#' @description Simulate GEDI full waveform data from Airborne Laser Scanning (ALS) 3D point cloud
#'
#' Input and output filenames, and formats
#' @param input character vector. lasfile input filename
# inList should be parsed from vector input
# @param inList list. input file list (ASCII file) for multiple files
#' @param output character. output filename

# Ground should always be true
#' @param ground record separate ground and canopy waveforms, default TRUE (shouldn't change).

# Always HDF
# @param hdf write output as HDF5. Best with gridded or list of coords
#' @param ascii write output as ASCII. Good for quick tests, default FALSE
#' @param waveID id. supply a waveID to pass to the output (only for single footprints)
#' Single footprint, list of footprints, or grid of footprints
#' @param coords lon lat numeric vector. footprint coordinate in same system as lasfile
#' @param listCoord name. Text file with list of coordinates. Pattern: X Y `[waveID]` `[geoCoordsX]` `[geoCoordsY]`. `[]` are optional, separated by spaces.
#' @param gridBound minX maxX minY maxY numeric vector. make a grid of waveforms in this box
#' @param gridStep res. grid step size
#' Lidar characteristics. Defaults are expected GEDI values.
#' @param pSigma pSigmasig. set Gaussian pulse width as 1 sigma
#' @param pFWHM fhwm. set Gaussian pulse width as FWHM in ns
#' @param readPulse file. read pulse shape and width from a file instead of making Gaussian
#' @param fSigma sig. set footprint width
#' @param wavefront file. read wavefront shape from file instead of setting Gaussian. Note that footprint width is still set by fSigma
#' @param res res. range resolution of waveform digitisation to output, in units of ALS data

# Not LVIS
# @param LVIS use LVIS pulse length, sigma=6.25m
#' @param topHat use a top hat wavefront
#' @param sideLobe use side lobes
#' @param lobeAng ang. lobe axis azimuth
#' Input data quality filters
#' @param checkCover check that the footprint is covered by ALS data. Do not output if not
#' @param maxScanAng ang. maximum scan angle, degrees
#' @param decimate x. probability of accepting an ALS beam
#' Computational speed options
#' @param pBuff s. point reading buffer size in Gbytes
#' @param maxBins for HDF5, limiting number of bins to save trimming.
#' @param countOnly only use count method
#' @param pulseAfter apply the pulse smoothing after binning for computational speed, at the risk of aliasing (default)
#' @param pulseBefore apply the pulse smoothing before binning to avoid the risk of aliasing, at the expense of computational speed
#' @param noNorm don't normalise for ALS density

#' Octree
#' @param noOctree do not use an octree
#' @param octLevels n. number of octree levels to use
#' @param nOctPix n. number of octree pixels along a side for the top level

#' Using full waveform input data (not tested)
# Not supported yet
# @param decon deconvolve
# @param indDecon deconvolve individual beams
# @param readWave read full waveform where available
# Miscellaneous

# R user should never only list files
# @param listFiles list files. Do not read them
#' @param keepOld do not overwrite old files, if they exist
#' @param useShadow account for shadowing in discrete return data through voxelization
#' @param polyGround find mean ground elevation and slope through fitting a polynomial

# nnGround is not working yet
# @param nnGround find mean ground elevation and slope through nearest neighbour
# @param seed n integer. random number seed
#'
#' #'
#' @return Returns an S4 object of class [`hdf5r::H5File-class`]
#' containing the simulated GEDI full-waveform.
#'
#' @seealso
#' i) Hancock, S., Armston, J., Hofton, M., Sun, X., Tang, H., Duncanson, L.I., Kellner,
#' J.R. and Dubayah, R., 2019. The GEDI simulator: A large-footprint waveform lidar simulator
#' for calibration and validation of spaceborne missions. Earth and Space Science.
#' \doi{10.1029/2018EA000506}
#'
#' ii) gediSimulator: \url{https://bitbucket.org/StevenHancock/gedisimulator/src/master/}
#'
#' @examples
#' \dontshow{
#' rm(list = ls())
#' }
#' libsAvailable <- require(lidR) && require(plot3D)
#' if (libsAvailable) {
#'     outdir <- tempdir()
#'
#'     # specify the path to ALS data (zip)
#'     alsfile_Amazon_zip <- system.file("extdata", "Amazon.zip", package = "rGEDIsimulator")
#'     alsfile_Savanna_zip <- system.file("extdata", "Savanna.zip", package = "rGEDIsimulator")
#'
#'     # Unzipping ALS data
#'     alsfile_Amazon_filepath <- unzip(alsfile_Amazon_zip, exdir = outdir)
#'     alsfile_Savanna_filepath <- unzip(alsfile_Savanna_zip, exdir = outdir)
#'
#'     # Reading and plot ALS file (las file)
#'     als_Amazon <- readLAS(alsfile_Amazon_filepath)
#'     als_Savanna <- readLAS(alsfile_Savanna_filepath)
#'
#'     # Extracting plot center geolocations
#'     xcenter_Amazon <- mean(st_bbox(als_Amazon)[c(1, 3)])
#'     ycenter_Amazon <- mean(st_bbox(als_Amazon)[c(2, 4)])
#'     xcenter_Savanna <- mean(st_bbox(als_Savanna)[c(1, 3)])
#'     ycenter_Savanna <- mean(st_bbox(als_Savanna)[c(2, 4)])
#'
#'     # Simulating GEDI full waveform
#'     wf_Amazon <- gediWFSimulator(
#'         input = alsfile_Amazon_filepath,
#'         output = file.path(outdir, "gediWF_amazon_simulation.h5"),
#'         coords = c(xcenter_Amazon, ycenter_Amazon)
#'     )
#'
#'     wf_Savanna <- gediWFSimulator(
#'         input = alsfile_Savanna_filepath,
#'         output = file.path(outdir, "gediWF_Savanna_simulation.h5"),
#'         coords = c(xcenter_Savanna, ycenter_Savanna)
#'     )
#'
#'     # Plot ALS and GEDI simulated full waveform
#'
#'     oldpar <- par()
#'     par(mfrow = c(2, 2), mar = c(4, 4, 0, 0), oma = c(0, 0, 1, 1), cex.axis = 1.2)
#'     scatter3D(als_Amazon@data$X, als_Amazon@data$Y, als_Amazon@data$Z,
#'         pch = 16, colkey = FALSE, main = "", cex = 0.5, bty = "u", 
#'         col.panel = "gray90", phi = 30, alpha = 1, theta = 45, col.grid = "gray50",
#'         xlab = "UTM Easting (m)", ylab = "UTM Northing (m)", zlab = "Elevation (m)"
#'     )
#'
#'     # Simulated waveforms shot_number is incremental beggining from 0
#'     shot_number <- 0
#'     simulated_waveform_amazon <- getLevel1BWF(wf_Amazon, shot_number)
#'     plot(simulated_waveform_amazon,
#'         relative = TRUE, polygon = TRUE, type = "l", lwd = 2, col = "forestgreen",
#'         xlab = "", ylab = "Elevation (m)", ylim = c(90, 140)
#'     )
#'     grid()
#'     scatter3D(als_Savanna@data$X, als_Savanna@data$Y, als_Savanna@data$Z,
#'         pch = 16, colkey = FALSE, main = "", cex = 0.5, bty = "u",
#'         col.panel = "gray90", phi = 30, alpha = 1, theta = 45, col.grid = "gray50",
#'         xlab = "UTM Easting (m)", ylab = "UTM Northing (m)", zlab = "Elevation (m)"
#'     )
#'
#'     shot_number <- 0
#'     simulated_waveform_savanna <- getLevel1BWF(wf_Savanna, shot_number)
#'     plot(simulated_waveform_savanna,
#'         relative = TRUE, polygon = TRUE, type = "l", lwd = 2, col = "green",
#'         xlab = "Waveform Amplitude (%)", ylab = "Elevation (m)", ylim = c(815, 835)
#'     )
#'     grid()
#'
#'     par(oldpar)
#'
#'     close(wf_Amazon)
#'     close(wf_Savanna)
#' }
#'
#' @import rGEDI hdf5r fs
#' @useDynLib rGEDIsimulator
#' @export
gediWFSimulator <- function(
    input,
    output,
    ground = TRUE,
    ascii = FALSE,
    waveID = NULL,
    coords = NULL,
    listCoord = NULL,
    gridBound = NULL,
    gridStep = 30.0,
    pSigma = -1.0,
    pFWHM = 15.0,
    readPulse = NULL,
    fSigma = 5.5,
    wavefront = NULL,
    res = 0.15,
    topHat = FALSE,
    sideLobe = FALSE,
    lobeAng = 0.0,
    checkCover = FALSE,
    maxScanAng = 1000000.0,
    decimate = 1.0,
    pBuff = as.integer(200000000),
    maxBins = as.integer(1024),
    countOnly = FALSE,
    pulseAfter = FALSE,
    pulseBefore = TRUE,
    noNorm = FALSE,
    noOctree = FALSE,
    octLevels = as.integer(0),
    nOctPix = as.integer(40),
    keepOld = FALSE,
    useShadow = FALSE,
    polyGround = FALSE) {
    # Set parameters that shouldn't be changed for GEDI
    l1b <- TRUE
    hdf <- FALSE
    LVIS <- FALSE
    if (ascii == TRUE) {
        hdf <- FALSE
        l1b <- FALSE
    }

    decon <- FALSE
    indDecon <- FALSE
    readWave <- FALSE

    listFiles <- FALSE
    nnGround <- FALSE

    # Check values
    stopifnotMessage(
        "input file(s) do not exist!" = all(file.exists(input)),
        "input is not a LAS file!" = all(fs::path_ext(input) == "las"),
        "output path is not valid!" = dir.exists(fs::path_dir(output)),
        "ascii should is not valid!" = checkLogical(ascii),
        "waveID should only work along with coords!" = is.null(waveID) || length(coords) == 2, # If waveID should only work along with coords
        "coords is invalid!" = checkNumericLength(coords, 2),
        "listCoord is invalid!" = is.null(listCoord) || file.exists(listCoord),
        "gridBound is invalid!" = checkNumericLength(gridBound, 4),
        "gridStep is invalid!" = checkNumeric(gridStep),
        "pFWHM is invalid!" = checkNumeric(pFWHM),
        "readPulse is invalid!" = checkFilepath(readPulse, newFile = FALSE, optional = TRUE),
        "fSigma is invalid!" = checkNumeric(fSigma),
        "wavefront is invalid!" = checkFilepath(wavefront, newFile = FALSE, optional = TRUE),
        "res is invalid!" = checkNumeric(res),
        "topHat is invalid!" = checkLogical(topHat),
        "sideLobe is invalid!" = checkLogical(sideLobe),
        "lobeAng is invalid!" = checkNumeric(lobeAng),
        "checkCover is invalid!" = checkLogical(checkCover),
        "maxScanAng is invalid!" = checkNumeric(maxScanAng),
        "decimate is invalid!" = checkNumeric(decimate),
        "pBuff is invalid!" = checkNumeric(pBuff),
        "maxBins is invalid!" = checkInteger(maxBins),
        "countOnly is invalid!" = checkLogical(countOnly),
        "pulseAfter is invalid!" = checkLogical(pulseAfter),
        "pulseBefore is invalid!" = checkLogical(pulseBefore),
        "noNorm is invalid!" = checkLogical(noNorm),
        "noOctree is invalid!" = checkLogical(noOctree),
        "octLevels is invalid!" = checkInteger(octLevels),
        "nOctPix is invalid!" = checkInteger(nOctPix),
        "keepOld is invalid!" = checkLogical(keepOld),
        "useShadow is invalid!" = checkLogical(useShadow),
        "polyGround is invalid!" = checkLogical(polyGround)
    )

    if (is.null(coords) && is.null(listCoord) && is.null(gridBound)) {
        stop("Coordinates for the waveforms should be provided!\nTIP: Use coords, listCoord or gridBound.")
    }

    inputInList <- inputOrInList(input)
    if (fs::path_ext(output) != "h5") {
        output <- fs::path_ext_set(output, ".h5")
    }

    .Call(
        "C_gediSimulator",
        inputInList[[1]],
        output,
        inputInList[[2]],
        ground,
        hdf,
        ascii,
        l1b,
        waveID,
        coords,
        listCoord,
        gridBound,
        gridStep,
        pSigma,
        pFWHM,
        readPulse,
        fSigma,
        wavefront,
        res,
        LVIS,
        topHat,
        sideLobe,
        lobeAng,
        checkCover,
        maxScanAng,
        decimate,
        as.integer(pBuff),
        as.integer(maxBins),
        countOnly,
        pulseAfter,
        pulseBefore,
        noNorm,
        noOctree,
        as.integer(octLevels),
        as.integer(nOctPix),
        decon,
        indDecon,
        readWave,
        listFiles,
        keepOld,
        useShadow,
        polyGround,
        nnGround,
        NULL
    )

    unloadLibrary()

    cleanInList(inputInList)
    if (ascii == TRUE) {
        result <- utils::read.table(output)
        metadata <- strsplit(readLines(output)[2:7], split = " ")
        final <- list()
        if (ncol(result) == 10) {
            colnames(result) <- c("elevation", "discrete intensity", "discrete count", "discrete fraction", "ALS pulse", "ALS and GEDI pulse", "ind decon", "ind decon GEDI", "decon GEDI", "ind decon")
        } else if (ncol(result) == 16) {
            colnames(result) <- c("elevation", "discrete intensity", "int canopy", "int ground", "discrete count", "count canopy", "count ground", "discrete fraction", "fraction canopy", "fraction ground", "ALS pulse", "ALS and GEDI pulse", "ind decon", "ind decon GEDI", "decon GEDI", "ind decon")
        }
        temp <- as.numeric(metadata[[1]][c(3, 5, 7, 9)])
        names(temp) <- c("fSigma", "pSigma", "res", "sideLobes")
        final <- as.list(temp)
        print(metadata)
        final$"coord" <- as.numeric(metadata[[2]][c(3, 4)])
        final$density_point <- as.numeric(metadata[[3]][4])
        final$beam <- as.numeric(metadata[[3]][6])
        final$meanScanAng <- as.numeric(metadata[[4]][3])
        line <- 5
        if (is.null(waveID) == FALSE) {
            final$waveID <- waveID
            line <- line + 1
        }
        final$ground <- as.numeric(metadata[[line]][3])
        final$groundAlg <- ifelse(polyGround, "poly", "simple")
        final$slope <- as.numeric(metadata[[line]][4])
        final$waveform <- result
        return(final)
    }
    result <- tryCatch(hdf5r::H5File$new(output, "r"), error = function(e) stop("The output file was not created\nSomething went wrong!"))

    result <- new("gedi.level1b", h5 = result)

    return(result)
}
