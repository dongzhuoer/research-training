% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{ref_species}
\alias{ref_species}
\title{return ref species for a species}
\usage{
ref_species(name)
}
\arguments{
\item{name}{character. see \link{extract_genus}}
}
\value{
character. ref_species paratmeter for HaMStR
}
\description{
return ref species for a species
}
\examples{
{
    aves::ref_species(aves::transcript$name);
    aves::ref_species(aves::outgroup$name);
}
}
