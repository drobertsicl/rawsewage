% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{timeTest}
\alias{timeTest}
\title{Estimate potential image dimensions by scan time disjunction}
\usage{
timeTest(input_vec, test_vec)
}
\arguments{
\item{input_vec}{A vector of a subset of linearly increasing scan times}

\item{test_vec}{The full length vector of scan times for the file in question}
}
\value{
Vector of scan time disjunction indixes
}
\description{
Searches for sudden jumps in scan acquisition time to calculate raster length
}
\details{
timeTest simply looks at the difference between sequential scan times in a given selection of scans vs the difference over the entire course of the run. Where acquisition stops for the stage to rasterise, this results in a sudden increase in scan time at the rasterisation position, which should be greater than the standard deviation.
}
\examples{
timeTest(scan_meta$scan_time[1:50], scan_meta$scan_time)
}
