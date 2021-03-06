\name{ICIResponse}
\alias{ICIResponse}
\title{
Markers and response/survival to ICI therapy
}
\description{
To explore the association between a marker and ICI therapy in datasets.
}
\usage{
ICIResponse(marker, Gide = Gide, Riaz = Riaz, Kim = Kim, Cho = Cho,
 Jung = Jung, background = background)
}
\arguments{
  \item{marker}{
if marker was a vector, we calculated the average of the marker genes;if marker was a list of signature genes, we use singscore method calculated the signature score. Your marker should contain a list of Upregulated and Downregulated (if necessary) gene-set names. 
}
  \item{Gide}{
expression matrix and clinical information from Gide et al, see Gide for more details.
}
  \item{Riaz}{
expression matrix and clinical information from Riaz et al, see Riaz for more details.
}
  \item{Kim}{
expression matrix and clinical information from Kim et al, see Kim for more details.
}
  \item{Cho}{
expression matrix and clinical information from Cho et al, see Cho for more details.
}
  \item{Jung}{
expression matrix and clinical information from Jung et al, see Jung for more details.
}
  \item{background}{
background gene sets for getSigScore, see background for more detail.
}
}
\details{
We return a list. It contains the relationship to ICI therapy response (a dataframe contains mean scores in response/non-response group, pvalues in Wilcox comparison and AUCs) and cox proportional model (continous variable).
}
\examples{
data(PD1hiCD8T_Sig)
data(Gide)
data(Riaz)
data(Kim)
data(Jung)
data(Cho)
data(background)
ICIResponseRes <- ICIResponse(list(PD1hiCD8T_Sig$Gene.Name), Gide = Gide, Riaz = Riaz,
 Kim = Kim,  Cho = Cho, Jung = Jung, background = background)
}
