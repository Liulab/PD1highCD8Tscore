\name{getSigScore}
\alias{getSigScore}
\title{Calculate PD-1hi CD8+ T cell signature score}
\description{
You can calculate PD-1hi CD8+ T cell signature score in data mentioned in our paper or your own sequencing data. Or your can calculate your own signature. 
}
\usage{
getSigScore(tpm, signature_genes, background = background)
}
\arguments{
  \item{tpm}{A gene expression matrix contains each gene (in rows) and sample (in columns).}
  \item{signature_genes}{A list of character vectors for gene sets, that can consist of gene IDs.}
  \item{background}{background gene sets for getSigScore, see background for more detail.}
}
\examples{
data(Cho)
data(background)
data(PD1hiCD8T_Sig)
getSigScoreRes = getSigScore(Cho$"RNA",list(PD1hiCD8T_Sig$Gene.Name),background = background)
}
