% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createMassBank.R
\name{gatherCCTE}
\alias{gatherCCTE}
\title{Retrieve supplemental annotation data from US EPA}
\usage{
gatherCCTE(key, api_key)
}
\arguments{
\item{key}{An Inchi-Key or other chemical identifier (e.g. Chemical name, DTXSID, CASRN, InChIKey, DTXCID)}

\item{api_key}{An US EPA CCTE API key (personal or application)}
}
\value{
Returns a list with 5 slots:
\code{dtxsid} The US EPA chemical dashboard substance id
\code{dtxcid} The US EPA chemical dashboard chemical id
\code{preferredName} The US EPA chemical dashboard preferred name
\code{casrn} The latest CAS registration number
\code{smiles} The SMILES annotation of the structure
}
\description{
Retrieves annotation data for a compound from the internet service US EPA CCTE 
based on the inchikey generated by babel or Cactus
}
\details{
The data retrieved is the US EPA DTXSID, the US EPA chemical dashboard
substance ID, the CAS-RN, the DTX preferred name, and the DTXCID (chemical ID).
}
\examples{

# Gather data for compound ID 131
\dontrun{gatherCCTE("QEIXBXXKTUNWDK-UHFFFAOYSA-N", api_key = NA)}

}
\references{
CCTE REST:
\url{https://api-ccte.epa.gov/docs/}
}
\seealso{
\code{\link{mbWorkflow}}
}
\author{
Tobias Schulze
}
