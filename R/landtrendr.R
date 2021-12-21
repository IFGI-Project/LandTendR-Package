## Creating Global Environment for Initial test
pkg.globals <- new.env()

## Initial message for Users
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("LandTrendr in  R.")
  packageStartupMessage("------------------------------------------")
  packageStartupMessage("This package execute LandtrendR algorithm based on GEE library")
  packageStartupMessage("ee_install will be run and a initialization of GEE is required")
}

## Global variables setup
##
## @description  In this section the global parameters for initial config are setup
##
## @param aoi_ Initial Area Of Interest
## @param distDir Parameter of Direction
## @param global_sensor parameter of the sensor to use (e.g L08)
## @param global_median_calcDif mean calculated based on the Collection
## @param indexNameFTV Initialize the FTV param
## @param runParam Initial config for segmentation process
## @param dummyCollection Initialize the a collection of zeros

rgee::ee_Initialize()
pkg.globals$time <- reticulate::import("time")
pkg.globals$datetime <- reticulate::import("datetime")
pkg.globals$filename <- system.file("shape/nc.shp", package="sf")
pkg.globals$nc <- sf::st_read(pkg.globals$filename)
pkg.globals$aoi <- sf::st_geometry(pkg.globals$nc[c(12),])
pkg.globals$aoi_ <- rgee::sf_as_ee(pkg.globals$aoi)

pkg.globals$distDir = -1
pkg.globals$global_sensor = NULL
pkg.globals$global_median_calcDif = NULL
pkg.globals$indexNameFTV = NULL
pkg.globals$runParam = list(
    maxSegments = 6,
    spikeThreshold = 0.9,
    vertexCountOvershoot = 3,
    preventOneYearRecovery =  TRUE,
    recoveryThreshold = 0.25,
    pvalThreshold = 0.05,
    bestModelProportion = 0.75,
    minObservationsNeeded= 6)

pkg.globals$dummyCollection = rgee::ee$ImageCollection(c(rgee::ee$Image(c(0,0,0,0,0,0))$mask(rgee::ee$Image(0))))

######################## CLASS METHODS ########################

#' Set a new AOI
#'
#' @description This function allows to set up a new Area of Interest
#'
#' @param aoi_ new Area of Interest, should be a sf geometry
#' @export
#'
setAOI <- function(aoi_){
  pkg.globals$aoi_ <- rgee::sf_as_ee(aoi_)
}

#' LandTrendR running Parameter
#'
#' @description The function will change the initial parameters to process the algorithm
#'
#' @param param List of parameter
#' @examples
#' list( maxSegments = 6,
#' spikeThreshold = 0.9,
#' vertexCountOvershoot = 3,
#' preventOneYearRecovery =  TRUE,
#' recoveryThreshold = 0.25,
#' pvalThreshold = 0.05,
#' bestModelProportion = 0.75,
#' minObservationsNeeded= 6)
#' @export
#'
setParameters <-  function(param)
{
  pkg.globals$runParam <- param
}

#' Harmonize L8 images
#'
#' @description This function harmonized the L8 image as the sensor has different specification as L5 and L7
#'
#' @param oli L8 Image in ee.Image format.
#' @return oli_harmonized ee.Image object
#'
harmonizationRoy <- function(oli){
  slopes = rgee::ee$Image$constant(c(0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949))
  itcp <- rgee::ee$Image$constant(c(-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029))
  oli_harmonized <- oli$select(c('B2','B3','B4','B5','B6','B7'),c('B1', 'B2', 'B3', 'B4', 'B5', 'B7'))$resample('bicubic')$subtract(itcp$multiply(10000))$divide(slopes)$set('system:time_start', oli$get('system:time_start'))$toShort()
  return(oli_harmonized)
}

#' Dynamic Sensor variable
#'
#' @description This function change dynamically all the sensor variables during the process.
#'
#' @param sensor string variable that refers to the sensor name in GEE syntaxis (i.e L08,L07,L05)
#' @export
set_global_sensor <- function(sensor){
  pkg.globals$global_sensor <- sensor
}

#' Image preparation function
#'
#' @description This function mask the images in the Image collection built with the LtRun input varibles.
#'
#' @param img ee.Image variable that refers to images in the preselected image collection based on input parameter od LtRun
#'
prepImages <- function(img){
  dat <- rgee::ee$Image(
    rgee::ee$Algorithms$If(
      pkg.globals$global_sensor == 'LC08',
      harmonizationRoy(img$unmask()),
      img$select(c('B1', 'B2', 'B3', 'B4', 'B5', 'B7'))$unmask()$resample('bicubic')$set('system:time_start', img$get('system:time_start'))
    )
  )
  qa <- img$select('pixel_qa')
  mask <- qa$bitwiseAnd(8)$eq(0)$And(qa$bitwiseAnd(16)$eq(0))$And(qa$bitwiseAnd(32)$eq(0))
  dat_masked <- dat$mask(mask)
  return(dat_masked)
}

#' SR Collection
#'
#' @description This function search for the images and create a ee.ImageCollection for each sentinel mission L05,L07,L08. This function also applies the Image preparation function
#'
#' @param year this parameter is automatically and refers to each year from Start Year and End Year
#' @param startDay this parameter refers to month and day of the starting year
#' @param endDay this parameter refers to month and day of the ending year
#' @param sensor this parameter refers sensor in this case either L08, L07 or L05
#' @param aoi this parameter refers Area Of Interest and constrain the retrieve images to the AOI.
#'
getSRcollection <- function(year, startDay, endDay, sensor, aoi){
  set_global_sensor(sensor)
  startDate <- paste(year,'-',startDay, sep = "")
  endDate <- paste(year,'-',endDay, sep = "")
  satellite <- paste("LANDSAT/",sensor,"/C01/T1_SR", sep="")
  srCollection <- rgee::ee$ImageCollection(satellite)$filterBounds(aoi)$filterDate(startDate, endDate)
  srCollection_prep <- srCollection$map(prepImages)
  return(srCollection_prep)
}

#' SR Collection merge
#'
#' @description This function merge the outcome Image collection from each Landsat mission.
#'
#' @param year this parameter is automatically and refers to each year from Start Year and End Year
#' @param startDay this parameter refers to month and day of the starting year
#' @param endDay this parameter refers to month and day of the ending year
#' @param aoi this parameter refers Area Of Interest and constrain the retrieve images to the AOI.
#'
getCombinedSRcollection <- function(year, startDay, endDay, aoi){
  lt5 = getSRcollection(year, startDay, endDay, 'LT05', aoi)
  le7 = getSRcollection(year, startDay, endDay, 'LE07', aoi)
  lc8 = getSRcollection(year, startDay, endDay, 'LC08', aoi)
  return (rgee::ee$ImageCollection(lt5$merge(le7)$merge(lc8)))
}

#' Global Median Difference
#'
#' @description This function change dynamically the global median different variable.
#'
#' @param finalCollection this parameter is automatically and refers to each year image in the construction of mosaics
#'
set_global_median_calcDif <- function(finalCollection){
  pkg.globals$global_median_calcDif <<- finalCollection$median()
}

#' Difference calculation from Medoid function
#'
#' @description This function change dynamically and calculates de Difference image from Medoid fuction
#'
#' @param img ee.Image
#' @return  reduced Image ee_Image
calcDifFromMed <- function(img){
  diff <- rgee::ee$Image(img)$subtract(pkg.globals$global_median_calcDif)$pow(rgee::ee$Image$constant(2))
  return(diff$reduce('sum')$addBands(img))
}

#' Medoid Function
#'
#' @description This function creates an Image collection with 6 bands
#'
#' @param inCollection ee.ImageCollection from the range of dates provided
#' @param dummyCollection empty ee.ImageCollection filled with zeros
#' @return Image Collection
#'
medoidMosaic <- function(inCollection, dummyCollection){
  imageCount <- inCollection$toList(1)$length()
  finalCollection <- rgee::ee$ImageCollection(rgee::ee$Algorithms$If(imageCount$gt(0), inCollection, dummyCollection))
  set_global_median_calcDif(finalCollection)
  difFromMedian <- finalCollection$map(calcDifFromMed)
  return(rgee::ee$ImageCollection(difFromMedian)$reduce(rgee::ee$Reducer$min(7))$select(c(1,2,3,4,5,6), c('B1','B2','B3','B4','B5','B7')))
}

#' Bluid Mosaic Function
#'
#' @description This function creates a mosaic from ee.Image
#' @param year this parameter is automatically and refers to each year from Start Year and End Year
#' @param startDay String variable - date tie to the Starting year (Default: '06-01')
#' @param endDay String variable - date tie to the Ending year (Default: '09-30')
#' @param aoi this parameter refers Area Of Interest and constrain the retrieve images to the AOI.
#' @param dummyCollection empty ee.ImageCollection filled with zeros
#' @return ee.Image
#'
buildMosaic <- function(year, startDay, endDay, aoi, dummyCollection){
  collection <- getCombinedSRcollection(year, startDay, endDay, aoi)
  time_Start <- pkg.globals$datetime$date(as.integer(year),as.integer(8),as.integer(1))
  img <- medoidMosaic(collection, dummyCollection)$set('system:time_start', as.numeric(as.POSIXct(time_Start)))
  return(rgee::ee$Image(img))
}

#' Bluid Mosaic Collection Function
#'
#' @description This function creates an Image collection with 6 bands
#'
#' @param startYear Numeric variable - Starting year of th analysis (Default: 2010)
#' @param endYear Numeric variable - Ending year of th analysis (Default: 2017)
#' @param startDay String variable - date tie to the Starting year (Default: '06-01')
#' @param endDay String variable - date tie to the Ending year (Default: '09-30')
#' @param aoi this parameter refers Area Of Interest and constrain the retrieve images to the AOI.
#' @param dummyCollection empty ee.ImageCollection filled with zeros
#'
#' @return ee.Image

buildMosaicCollection <- function(startYear, endYear, startDay, endDay, aoi, dummyCollection){
  imgs = list()
  for (year in seq(startYear,endYear+1)){
    tmp <- buildMosaic(year, startDay, endDay, aoi, dummyCollection)
    time_Start <- pkg.globals$datetime$date(as.integer(year),as.integer(8),as.integer(1))
    imgs <-c(imgs,tmp$set('system:time_start', as.numeric(as.POSIXct(time_Start))))
  }
  return(rgee::ee$ImageCollection(imgs))
}

#' segIndex Function
#'
#' @description This function calculates and add to the Image the calculation of NBR index
#'
#' @param img ee.Image
#' @return ee.Image
#'
segIndex <- function(img){
  index <- img$normalizedDifference(c('B4', 'B7'))$multiply(1000)$select(0)$rename("NBR")$set('system:time_start', img$get('system:time_start'))
  return(index$toShort())
}

#' Invert index Function
#'
#' @description This function invert the direction of the magnitud to make it logic for the user.
#'
#' @param img ee.Image
#' @return ee.Image

invertIndex <- function(img){
  return(img$multiply(pkg.globals$distDir)$toShort()$set('system:time_start', img$get('system:time_start')))

}

#' Invert index Function
#'
#' @description This function invert the FTV.
#'
#' @param img ee.Image
#' @return ee.Image
#'
invertFTV <- function(img){
  return(img$addBands(img$select(0)$rename(pkg.globals$indexNameFTV)$multiply(pkg.globals$distDir))$toShort()$set('system:time_start', img$get('system:time_start')))
}

#################################################### EXPORT FUNCTION DEFINITION ########################################
#' LTvertStack Function
#'
#' @description This function extracts the segmentation vertex info
#'
#' @param LTresult LandtrendR result
#' @return ee.Image
#' @export

getLTvertStack <- function(LTresult){
  emptyArray <- list()
  vertLabels <- list()
  for(i in seq(1, pkg.globals$runParam['maxSegments'][[1]]+2)){
    label <- paste("vert_",as.character(i),sep = "")
    vertLabels <- c(vertLabels, label)
  }
  emptyArray <- c(emptyArray,0)

  zeros <- rgee::ee$Image(rgee::ee$Array(c(emptyArray,
                                           emptyArray,
                                           emptyArray)))

  lbls <- c(c('yrs_','src_','fit_'), c(vertLabels))

  vmask <- LTresult$arraySlice(0,3,4)

  ltVertStack <- LTresult$arrayMask(vmask)$arraySlice(0, 0, 3)$addBands(zeros)$toArray(1)

  return(ltVertStack)
}

#' Run LandTrendr Algorithm
#'
#' This function built a ImageCollection suitable for LT analysis, using Harmonization and pre-processes before running the Landtrendr Algorithms
#' The function starts building the image and an additional NBR index that helps the algorithm to calculate the magnitude of change.
#' This function is based on rgee package, thus the outcome of this function is a ee.Image object.
#'
#' @param startYear Numeric variable - Starting year of th analysis (Default: 2010)
#' @param endYear Numeric variable - Ending year of th analysis (Default: 2017)
#' @param startDay String variable - date tie to the Starting year (Default: '06-01')
#' @param endDay String variable - date tie to the Ending year (Default: '09-30')
#' @param aoi Polygon sf Geometry - Area of Interest (Default: US County from sf dataset)
#' @param runParams Parameters for Landtrendr Analysis - By default the  package take the following values: maxSegments = 6, spikeThreshold = 0.9, vertexCountOvershoot = 3, preventOneYearRecovery =  TRUE,recoveryThreshold = 0.25, pvalThreshold = 0.05, bestModelProportion = 0.75, minObservationsNeeded= 6
#' @return lt ee.Image object
#' @export
#'
LtRun <- function(startYear = 2010, endYear = 2017, startDay = '06-01', endDay = '09-30', aoi= pkg.globals$aoi_, runParams = pkg.globals$runParam){
  rgee::ee_Initialize()
  # build annual image collection and run LandTrendr
  annualSRcollection <- buildMosaicCollection(startYear, endYear, startDay, endDay, aoi, pkg.globals$dummyCollection)

  annualIndexCollection <- annualSRcollection$map(segIndex)$map(invertIndex)

  pkg.globals$indexNameFTV <<- paste("NBR","_FTV", sep = "")

  ltCollection = annualIndexCollection$map(invertFTV)

  # add the collection to the LandTrendr parameters and run LT-GEE
  runParams <- c(runParams,timeSeries=ltCollection)

  lt <- rgee::ee$Algorithms$TemporalSegmentation$LandTrendr(ltCollection,
                                                            runParams['maxSegments'][[1]],
                                                            runParams['spikeThreshold'][[1]],
                                                            runParams['vertexCountOvershoot'][[1]],
                                                            runParams['preventOneYearRecovery'][[1]],
                                                            runParams['recoveryThreshold'][[1]],
                                                            runParams['pvalThreshold'][[1]],
                                                            runParams['bestModelProportion'][[1]],
                                                            runParams['minObservationsNeeded'][[1]])
  return(lt)
}
