% Please edit documentation in R/dmlalg.R
\docType{package}
\name{dmlalg}
\alias{dmlalg}
\title{dmlalg: double machine learning algorithms}
\description{
The \code{dmlalg} package contains implementations of
double machine learning (DML) algorithms in \code{R}.
}
\section{Partially linear models with confounding variables}{
    Our goal is to perform inference for the linear parameter in partially
    linear models with confounding variables.
    The standard DML estimator of the linear parameter has a two-stage least
    squares interpretation, which can lead to a large variance and overwide
    confidence intervals.
    We apply regularization to reduce the variance of the estimator,
    which produces narrower confidence intervals that are approximately valid.
    Nuisance terms can be flexibly estimated with machine learning algorithms.
    \describe{
    \item{\code{\link{regsdml}}}{Estimates the
    linear parameter in a partially linear model with
    confounding variables with regularized and standard
    DML methods.}
    \item{\code{\link{summary.regsdml}}}{A \code{summary} method for objects
    fitted with \code{\link{regsdml}}.}
    \item{\code{\link{confint.regsdml}}}{A \code{confint} method for objects
    fitted with \code{\link{regsdml}}.}
    \item{\code{\link{coef.regsdml}}}{A \code{coef} method for objects
    fitted with \code{\link{regsdml}}.}
    \item{\code{\link{vcov.regsdml}}}{A \code{vcov} method for objects
    fitted with \code{\link{regsdml}}.}
    \item{\code{\link{print.regsdml}}}{A \code{print} method for objects
    fitted with \code{\link{regsdml}}.}
    }
}

\references{
C. Emmenegger and P. \enc{Bühlmann}{Buehlmann}.
Regularizing Double Machine Learning in Partially Linear Endogenous
Models, 2021. Preprint arXiv:2101.12525.
}

