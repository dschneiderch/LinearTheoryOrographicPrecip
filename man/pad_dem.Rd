% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pad_dem.R
\name{pad_dem}
\alias{pad_dem}
\title{Pad DEM with 0 to remove edge effects from FFT ----}
\usage{
pad_dem(h, pad = NULL)
}
\arguments{
\item{h}{matrix of topo}

\item{pad}{(optional) 0 for no padding. Leave blank for automatic padding.}
}
\description{
Pad DEM with 0 to remove edge effects from FFT ----
}
\details{
This will also remove any NA and convert to 0
}
