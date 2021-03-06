\name{getIPS}
\alias{getIPS}
\title{
Calculation of Immunophenoscore
}
\description{
Immunophenoscore (IPS) was calculated by a script supplied by The Cancer Immunome Atlas
}
\usage{
getIPS(gene_expression, IPSG = IPSG)
}
\arguments{
  \item{gene_expression}{ A gene expression matrix contains each gene (in rows) and sample (in columns).}
  \item{IPSG}{gene set in Immunophenoscore}
}
\references{
Charoentong P., Finotello F., Angelova M., et al. Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade. Cell Rep. Jan 3 2017;18(1):248-262.
}
\examples{
data(IPSG)
data(Cho)
getIPSRes <-getIPS(Cho$"RNA",IPSG = IPSG)
}
