defaultAOI <- function(){
  filename <- system.file("shape/nc.shp", package="sf")
  nc <- sf::st_read(filename)
  sel <- c(12)
  aoi_ <- sf::st_geometry(nc[sel,])
  aoi_ <- rgee::sf_as_ee(aoi_)
  return(aoi_)
}

globalParameters <-  function()
{
  me <- list(
    aoi_ = defaultAOI(),
    distDir = -1,
    global_sensor = NULL,
    global_median_calcDif = NULL,
    indexNameFTV = NULL,
    runParam = list(
      maxSegments = 6,
      spikeThreshold = 0.9,
      vertexCountOvershoot = 3,
      preventOneYearRecovery =  TRUE,
      recoveryThreshold = 0.25,
      pvalThreshold = 0.05,
      bestModelProportion = 0.75,
      minObservationsNeeded= 6),
    dummyCollection = rgee::ee$ImageCollection(c(rgee::ee$Image(c(0,0,0,0,0,0))$mask(rgee::ee$Image(0))))
  )
  ## Set the name for the class
  class(me) <- append(class(me),"globalParameters")
  return(me)
}

