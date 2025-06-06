% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{getMeta}
\alias{getMeta}
\title{Get metadata from a Waters .IDX metadata file}
\usage{
getMeta(idx_path, calib_coefs = NULL)
}
\arguments{
\item{idx_path}{A file path to a .IDX file}

\item{calib_coefs}{A vector of calibration coefficients computed by getCoefs}
}
\value{
Data frame
}
\description{
Automatically extracts vector of memory addresses indicating the start of each scan in an associated .DAT file
}
\examples{
getMeta("C:/My Folder/image.raw/_FUNC001.IDX")
}
