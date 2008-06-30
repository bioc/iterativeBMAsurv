\name{trainData}
\docType{data}
\alias{trainData}
\title{Sample Training Data for the Iterative BMA Algorithm for Survival Analysis}
\description{
This is an adapted diffuse large B-cell lymphoma (DLBCL) dataset from 
Rosenwald et al. (2002). This data matrix consists of the expression 
levels from 65 DLBCL samples (rows), and 100 top univariate genes 
(columns). This dataset is used as a sample training set in our examples.
}
\usage{data(trainData)}
\format{The data matrix is called \code{trainData}. Each entry
in the matrix represents the expression level of one gene from 
a DLBCL sample.
}
\details{
We started with the full expression data from Rosenwald et al. (2002),
which is available along with corresponding patient information at their 
supplemental website \url{http://llmpp.nih.gov/DLBCL/}. We selected
a subset of the 160 training samples, and then performed Cox Proportional 
Hazards Regression to obtain the 100 genes with the highest log likelihood. 
The filtered dataset consists of 65 training samples, and it is available at
our website \url{http://expression.washington.edu/ibmasurv/protected}.
} 
\source{
Full dataset: \url{http://llmpp.nih.gov/DLBCL/}.
}
\references{
Rosenwald, A., Wright, G., Wing, C., Connors, J., Campo, E. et al. (2002). 
The Use of Molecular Profiling to Predict Survival After Chemotherapy for Diffuse Large-B-Cell Lymphoma. 
The New England Journal of Medicine, 346(25), 1937-1947. 
}
\keyword{datasets}

