% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/landtrendr.R
\name{getSRcollection}
\alias{getSRcollection}
\title{SR Collection}
\usage{
getSRcollection(year, startDay, endDay, sensor, aoi)
}
\arguments{
\item{year}{this parameter is automatically and refers to each year from Start Year and End Year}

\item{startDay}{this parameter refers to month and day of the starting year}

\item{endDay}{this parameter refers to month and day of the ending year}

\item{sensor}{this parameter refers sensor in this case either L08, L07 or L05}

\item{aoi}{this parameter refers Area Of Interest and constrain the retrieve images to the AOI.}
}
\description{
This function search for the images and create a ee.ImageCollection for each sentinel mission L05,L07,L08. This function also applies the Image preparation function
}
