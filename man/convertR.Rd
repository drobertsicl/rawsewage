% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{convertR}
\alias{convertR}
\title{Convert Waters .raw file to .mzML via proteowizard}
\usage{
convertR(filename = NULL, outpath = NULL, msconvert = msconexe)
}
\arguments{
\item{filename}{Path to a .raw folder of imaging data}

\item{outpath}{Path for .mzML output to be written to}

\item{msconvert}{Path to msconvert.exe}
}
\value{
Returns nothing in the session, outputs a .mzML file somewhere on your hard drive, probably
}
\description{
Convenience function interfacing to proteowizard so the vignette can look relatively clean
}
\examples{
convertR(filename = "C:/My Folder/image.raw/", outpath = "C:/My Folder/", msconvert = "C:/Users/whoever/AppData/Local/Apps/ProteoWizard 3.0.24094.d2966db 64-bit/msconvert.exe")
}
