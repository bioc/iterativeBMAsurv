# Goal:      Compute the risk score of a given patient sample
# Arguments: newdata.vec = vector of test data (1 observation)
#            postprob.vec = postprob from bic.surv
#	     mle.mat = mle matrix from bic.surv (models = rows, variables = columns)
# Output:    A patient sample risk score

predictBicSurv <- function(newdata.vec, postprob.vec, mle.mat) {
       # First compute the risk score for each of the models
       # Risk score for each model = sum(coefficient in model k * corresponding expression level in newdata.vec)
       risk.score.vec <- apply (mle.mat, 1, function(x) sum(x * newdata.vec))
       
       # Patient sample risk score under BMA = sum (postprob for model k * risk score for model k)
       # over all the selected models
       retprob <- sum (postprob.vec * risk.score.vec)
       retprob
}
