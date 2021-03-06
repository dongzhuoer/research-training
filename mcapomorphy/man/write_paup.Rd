% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paup.R
\name{write_paup}
\alias{write_paup}
\title{perpare nexus file to run paup}
\usage{
write_paup(fasta_file, newick_file, outgroup, output_file)
}
\arguments{
\item{fasta_file}{string. sequence file in FASTA format. see
\code{\link[=fasta_to_nexus]{fasta_to_nexus()}} for requirement of the content}

\item{newick_file}{string. tree file in Newick format.}

\item{outgroup}{character. name of outgroups, \code{' '} would be replaced by
\code{'_'}}

\item{output_file}{string. path to the output NEXUS file.}
}
\value{
\code{NULL}
}
\description{
perpare nexus file to run paup
}
\details{
for \code{fasta_file}, \code{newick_file} and \code{outgroup} only
\code{[[:alnum:]_]} are allowed in taxon name (sequence header)
}
\examples{
temp_file <- tempfile()
write_paup(
    system.file('extdata', 'EOG090F05Z3.fasta', package = 'mcapomorphy'), 
    system.file('extdata', 'omics.tre', package = 'mcapomorphy'), 
    aves::outgroup$species, temp_file
)
\dontrun{
    system2('paup4a159', 'data-raw/EOG090F05Z3.nexus -n -u', T)
}

}
