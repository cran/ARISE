\name{ARISEFS}
\alias{ARISEFS}
\title{Adaptive Residual-Informed Stability Ensemble Feature Selection}
\usage{
ARISEFS(X, Y, topK = 10)
}
\arguments{
\item{X}{A real numeric or logical matrix, or a data frame with real numeric or logical vector columns. Rows are samples and columns are features. If supplied, feature names must be unique, non-missing, and non-empty.}
\item{Y}{An atomic response vector or factor with one non-missing class label for each row of \code{X}. At least two classes and at least three observations per class are required. Complex values and non-finite numeric labels are not supported.}
\item{topK}{A positive whole number giving the maximum number of top-ranked features to return. Defaults to 10.}
}
\description{
Selects the top-ranked features using the Adaptive Residual-Informed Stability Ensemble method.
}
\details{
ARISE combines seven relevance components with class-balanced subsample stability, within-class residual redundancy control, multiclass pair coverage, nested profile racing, and soft consensus. The feature path is ranked from strongest to weakest, and the requested leading feature indices are returned.

Missing and non-finite predictor values are replaced by column medians computed from the training data. Columns with no finite values or variance at most \code{1e-12} are removed. Preprocessing is fitted separately on each inner training split and on the complete input for final selection. An informative error is raised when no variable feature remains in a required training split.

Profile selection uses three stratified inner folds with fixed k-nearest neighbour, radial support vector machine, and random forest learners. The internal seed is fixed at 1; the caller's random number state is restored after internal randomized operations. Results are reproducible for the same input, R version, and dependency versions.

For an unbiased assessment of a prediction model, call this function separately on each outer training split. Select features before fitting the corresponding model and apply the returned column indices to its held-out data.
}
\value{
An unnamed integer vector containing unique column indices of the selected features in \code{X}, ordered from strongest to weakest. The length is at most \code{topK} and may be smaller when fewer usable features remain. No fitted classifier is returned.
}
\references{
Khan, Zardad, Amjad Ali, Naz Gul, and Saeed Aldahmani (2026).
ARISE: An adaptive residual-informed stability ensemble for feature selection
in small-sample biomedical omics. \emph{arXiv preprint arXiv:2608.14866}.
\doi{10.48550/arXiv.2608.14866}.
}
\examples{
X <- cbind(
  signal = c(-3, -2, -1, 1, 2, 3),
  noise = c(1, 3, 2, 3, 1, 2)
)
Y <- factor(rep(c("a", "b"), each = 3L))
selected <- ARISEFS(X, Y, topK = 1)
X[, selected, drop = FALSE]
}
\keyword{methods}
\keyword{multivariate}
