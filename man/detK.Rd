% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detK.R
\name{detK}
\alias{detK}
\title{determine the number of clusters}
\usage{
detK(D, epi = 1)
}
\arguments{
\item{D}{distance matrix}

\item{epi}{threshold for cutting off}
}
\value{
number of clusters
}
\description{
determine the number of clusters
}
\examples{
data(sim_dat)
data_counts = assays(sim_dat)$count
D_c = cal_D(data_counts,4)
detK(D_c)
}
