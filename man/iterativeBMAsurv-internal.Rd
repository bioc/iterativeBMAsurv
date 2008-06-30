\name{iterativeBMAsurv-internal}
\alias{iterativeBMAsurv-internal}
\alias{iterateBMAinit}
\alias{assignRiskGroup}
\alias{crossVal.init}
\alias{crossVal.run}
\alias{crossVal.fold}
\alias{crossVal.tabulate}
\alias{crossVal.final.calc}
\alias{imageplot.bma.mod}

\title{Internal functions for iterativeBMAsurv}
\description{Internal functions for \code{iterativeBMAsurv}, not meant to be called directly.}

\usage{
iterateBMAinit (x, maxNvar = 25)
assignRiskGroup (x, quantile.cutPoint)
crossVal.init (exset, noFolds, noRuns, nbest, maxNvar, p, cutPoint)
crossVal.run (exset, survTime, censor, diseaseType, nbest, maxNvar, p, cutPoint, rp.mat, overall.average.results, p.val.result.mat, chi.sq.result.mat, noFolds, cv, run, verbose)
crossVal.fold (cvpass, survtimepass, censorpass, diseaseType, nbest, maxNvar, p, cutPoint, lastSample, noFolds, single.run.result.mat, p.val.result.mat, chi.sq.result.mat, cv, fold, run, verbose)
crossVal.tabulate (results, single.run.result.mat)
crossVal.final.calc (overall.average.results, noRuns)
imageplot.bma.mod (bicreg.out, color = "default", ...)
}
\keyword{internal}

