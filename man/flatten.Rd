% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createMassBank.R
\name{flatten}
\alias{flatten}
\alias{readMbdata}
\title{Flatten, or re-read, MassBank header blocks}
\usage{
flatten(mbdata) 

readMbdata(row)
}
\arguments{
\item{mbdata}{A list of MassBank compound information sets as returned from
\code{\link{gatherData}}.}

\item{row}{One row of MassBank compound information retrieved from an
infolist.}
}
\value{
\code{flatten} returns a tibble (not a data frame or matrix) to be written to
CSV.

\code{readMbdata} returns a list of type \code{list(id= \var{compoundID},
..., 'ACCESSION' = '', 'RECORD_TITLE' = '', )} etc.
}
\description{
\code{flatten} converts a list of MassBank compound information sets (as
retrieved by \code{\link{gatherData}}) to a flat table, to be exported into
an \link[=loadInfolist]{infolist}. \code{readMbdata} reads a single record
from an infolist flat table back into a MassBank (half-)entry.
}
\details{
Neither the flattening system itself nor the implementation are particularly
fantastic, but since hand-checking of records is a necessary evil, there is
currently no alternative (short of coding a complete GUI for this and
working directly on the records.)
}
\examples{
\dontrun{
	# Collect some data to flatten
	ids <- c(40,50,60,70)
 data <- lapply(ids, gatherData)
 # Flatten the data trees to a table
 flat.table <- flatten(data)
 # reimport the table into a tree
 data.reimported <- apply(flat.table, 1, readMbdata)
}

}
\references{
MassBank record format:
\url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
}
\seealso{
\code{\link{gatherData}},\code{\link{loadInfolist}}
}
\author{
Michael Stravs
}
