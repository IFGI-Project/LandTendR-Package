library(rgee)
library(sf)
library(reticulate)

ee_Initialize()

# Geometry Parameter
filename <- system.file("shape/nc.shp", package="sf")
nc <- st_read(filename)
sel <- c(12)
aoi <- st_geometry(nc[sel,]) %>% sf_as_ee()
aoi

# Defining the Dates
startYear = 1985
endYear   = 2017
startDay  = '06-01'
endDay    = '09-30'

#Segmentation Index Function
segIndex <- function(img){
  index = img$normalizedDifference(c('B4', 'B7'))$multiply(1000)$select(c(0), C('NBR'))$set('system:time_start', img$get('system:time_start'))
  return(index$toShort())
}

# Global parameter definition

distDir <- -1

set_distDir <- function(dist){
    distDir <- dist
}

# Definition of runParameters
runParams = list(
  maxSegments = 6,
  spikeThreshold = 0.9,
  vertexCountOvershoot = 3,
  preventOneYearRecovery = 'true',
  recoveryThreshold = 0.25,
  pvalThreshold = 0.05,
  bestModelProportion = 0.75,
  minObservationsNeeded= 6)

# dummyCollection

dummyCollection = ee$ImageCollection(c(ee$Image(c(0,0,0,0,0,0))$mask(ee$Image(0))))

# Additional Methods Harmonization, SRCollection, CombinedSRcollection, medoidMosaic

harmonizationRoy <- function(oli){
   slopes = ee.Image.constant(c(0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949))
  itcp <- ee$Image$constant(c(-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029))
  return(oli$select(c('B2','B3','B4','B5','B6','B7'),c('B1', 'B2', 'B3', 'B4', 'B5', 'B7'))$resample('bicubic')$subtract(itcp$multiply(10000))$divide(slopes)$set('system:time_start', oli$get('system:time_start'))$toShort())
}

prepImages <- function(img){
  dat = ee$Image(
    ee$Algorithms$If(
      sensor == 'LC08',
      harmonizationRoy(img.unmask()),
      img$select(c('B1', 'B2', 'B3', 'B4', 'B5', 'B7'))
      $unmask()
      $resample('bicubic')
      $set('system:time_start', img$get('system:time_start'))
    )
  )
  qa = img$select('pixel_qa')
  mask = qa$bitwiseAnd(8)$eq(0)$And(qa$bitwiseAnd(16)$eq(0))$And(qa$bitwiseAnd(32)$eq(0))
  return(dat.mask(mask))
}
######################
getSRcollection <- function(year, startDay, endDay, sensor, aoi){
  srCollection = ee.ImageCollection('LANDSAT/'+ sensor + '/C01/T1_SR')$filterBounds(aoi)$filterDate(str(year)+'-'+startDay, str(year)+'-'+endDay)
  return(srCollection$map(prepImages))
}

getCombinedSRcollection <- function(year, startDay, endDay, aoi){
  lt5 = getSRcollection(year, startDay, endDay, 'LT05', aoi)
  le7 = getSRcollection(year, startDay, endDay, 'LE07', aoi)
  lc8 = getSRcollection(year, startDay, endDay, 'LC08', aoi)
  return (ee$ImageCollection(lt5$merge(le7)$merge(lc8)))
}
######################
calcDifFromMed <- function(img){
  diff = ee$Image(img)$subtract(median)$pow(ee.Image.constant(2))
  return(diff$reduce('sum')$addBands(img))
}

medoidMosaic <- function(inCollection, dummyCollection){
  imageCount = inCollection$toList(1)$length()
  finalCollection = ee.ImageCollection(ee.Algorithms.If(imageCount.gt(0), inCollection, dummyCollection))
  median = finalCollection.median()
  difFromMedian = finalCollection.map(calcDifFromMed)
  return(ee$ImageCollection(difFromMedian)$reduce(ee$Reducer$min(7))$select(c(1,2,3,4,5,6), c('B1','B2','B3','B4','B5','B7')))
}
#####################
buildMosaic <- function(year, startDay, endDay, aoi, dummyCollection){
  collection <- getCombinedSRcollection(year, startDay, endDay, aoi)
  img = medoidMosaic(collection, dummyCollection)$set('system:time_start', int(time.mktime(datetime.date(year,8,1)$timetuple())))
  return(ee$Image(img))
}

buildMosaicCollection <- function(startYear, endYear, startDay, endDay, aoi, dummyCollection){
  imgs = list()
  for (year in range(startYear,endYear+1)){
    tmp = buildMosaic(year, startDay, endDay, aoi, dummyCollection)
    imgs.append(tmp.set('system:time_start', int(time.mktime(datetime.date(year,8,1)$timetuple()))))
  }
  return(ee$ImageCollection(imgs))
}

getLTvertStack <- function(LTresult){
  emptyArray <- list()
  vertLabels <- list()
  for(i in range(1, runParams['maxSegments']+2)){
    vertLabels$append("vert_"+str(i))
  }
  emptyArray$append(0)

  zeros <- ee$Image(ee$Array(c(emptyArray,
                             emptyArray,
                             emptyArray)))

  lbls <- cbind(c('yrs_','src_','fit_'), c(vertLabels),)

  vmask <- LTresult$arraySlice(0,3,4)

  ltVertStack <- LTresult$arrayMask(vmask)$arraySlice(0, 0, 3)$addBands(zeros)$toArray(1)$arraySlice(1, 0, runParams['maxSegments']+1)$arrayFlatten(lbls, '')

  return(ltVertStack)
}

invertIndex <- function(img){
  return(img$multiply(distDir)$toShort()$set('system:time_start', img$get('system:time_start')))

}

invertFTV <- function(img){
  return(img.addBands(img$select(c(0),c(indexNameFTV))$multiply(distDir))$toShort()$set('system:time_start', img$get('system:time_start')))
}

##########################################################################

# build annual image collection and run LandTrendr
annualSRcollection <- buildMosaicCollection(startYear, endYear, startDay, endDay, aoi, dummyCollection)
annualIndexCollection <- annualSRcollection.map(segIndex).map(invertIndex)

indexNameFTV = ee$Image(annualIndexCollection$first())$bandNames()$getInfo()[0]+'_FTV'
ltCollection = annualIndexCollection$map(invertFTV)

# add the collection to the LandTrendr parameters and run LT-GEE
runParams <- c(runParams,timeSeries='ltCollection')
lt = ee$Algorithms$TemporalSegmentation$LandTrendr(runParams)
lt

ltVertStack <- getLTvertStack(lt$select("LandTrendr"))$toShort()

#Map$centerObject(eeObject = aoi)
#Map$addLayer(aoi)
