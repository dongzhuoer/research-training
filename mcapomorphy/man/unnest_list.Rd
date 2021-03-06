% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{unnest_list}
\alias{unnest_list}
\title{get the outermost list which contains more than one element}
\usage{
unnest_list(x)
}
\description{
\code{unnest_list(list(list(list(1))))} gives \code{list(1)}
}
\examples{
unnest_list(1)
unnest_list(list(1))
unnest_list(list(list(list(1))))
unnest_list(list(list(list(1), list(2))))
unnest_list(list(list(list(1)), list(2)))

}
