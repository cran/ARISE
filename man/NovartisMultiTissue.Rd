\name{NovartisMultiTissue}
\alias{NovartisMultiTissue}
\docType{data}
\title{Novartis Multitissue Gene Expression Dataset}
\description{
Gene expression data for 102 samples and 5,565 features from four tissue classes.
}
\usage{
data(NovartisMultiTissue)
}
\format{
A list with two components:
\describe{
\item{X}{A numeric matrix with 102 samples in rows and 5,565 features in columns.}
\item{Y}{A factor with 102 class labels. The levels \code{1}, \code{2}, \code{3}, and \code{4} contain 25, 26, 28, and 23 samples, respectively.}
}
}
\details{
The matrix contains processed expression values. Sample and feature identifiers are generated as \code{Sample1}, \code{Sample2}, \ldots{} and \code{Feature1}, \code{Feature2}, \ldots{}, respectively. The supplied file does not include mappings to tissue names, gene symbols, or probe identifiers. Two columns are constant and are removed during ARISE preprocessing.
}
\source{
Processed benchmark supplied with the ARISE source package. The associated original study is Su et al. (2002), cited below. No additional normalization is performed when loading the data.
}
\references{
Su, A. I., et al. (2002). Large-scale analysis of the human and mouse transcriptomes. \emph{Proceedings of the National Academy of Sciences of the United States of America}, 99(7), 4465--4470. \doi{10.1073/pnas.012025199}.
}
\examples{
data(NovartisMultiTissue)
X <- NovartisMultiTissue$X
Y <- NovartisMultiTissue$Y
}
\keyword{datasets}
