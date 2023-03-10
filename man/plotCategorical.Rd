% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotCategorical}
\alias{plotCategorical}
\title{Barplot wrapper}
\usage{
plotCategorical(
  ntinfo,
  field,
  ntID = NULL,
  na.rm = FALSE,
  main = NULL,
  cex = 0.5,
  file = NULL,
  width = 15,
  height = 15,
  bg = "white",
  units = "cm",
  res = 200
)
}
\arguments{
\item{ntinfo}{A data.frame with the input data. It should contain the
columns with the desired categorical data and a column labeled ntID.}

\item{field}{The column name with the desired data.}

\item{ntID}{A vector of integers with the desired nucleotides of
analysis. If NULL all the nucleotides in the data.frame will be used.}

\item{na.rm}{A logical to remove missing data.}

\item{main}{A string with the title of the plot.}

\item{cex}{To be passed to the par() function.}

\item{file}{A string with the name of the output file. If NULL, the
plot will be printed to screen.}

\item{width}{The width of the plot (passed to the png() function).}

\item{height}{The height of the plot (passed to the png() function).}

\item{bg}{The background color of the plot (passed to the png() function).}

\item{units}{The unit to measure height and width (passed to the png()
function).}

\item{res}{Resolution (passed to the png() function).}
}
\value{
A barplot with the categorical data of interest, which can be
directly saved  to a ".png" file.
}
\description{
Function to make more straigtforward the process of ploting a barplot for
categorical data.
}
\examples{
    ## To see all the types of trinucleotides in the dataset:
    ntinfo <- pipeNucData("1bau")
    plotCategorical(ntinfo=ntinfo, field="localenv")

}
\author{
Diego Gallego
}
