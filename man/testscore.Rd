\name{testscore}
\alias{testscore}
\title{Restrained adaptive test based on the standardized score process}
\description{
This function provides the statistic and the pvalue of the restrained adaptive test in Chauvel and OQuigley (2014).}
\usage{
testscore(formula, data, beta0=0, n_rep=10^6, digits=5)
}
\arguments{
  \item{formula}{
A formula object or character string with the time and censoring status separated by "+" on the left hand side and the covariates separated by "+" on the right. For instance, if the time name is "Time", the censoring status is "Status" and the covariates are "Cov1" and "Cov2", the formula is "Time+Status~Cov1+Cov2". No interaction can be provided.
}
  \item{data}{ 
  A data frame with the data. The censoring status should be 1 for failure and 0 for censoring. No missing data are accepted. } 
  \item{beta0}{
A vector of parameters to test in the null hypothesis H_0: beta = beta0. By default, beta0 = 0. Its length is the number of covariates. Each value corresponds to the regression coefficient for a covariate, in the same order as appearing in formula. 
}  
 \item{n_rep}{ 
  An integer for the number of simulations for the estimation of the p-value. It must be higher than 10^5.} 
 \item{digits}{ 
  An integer for the number of decimal places to be used in the results.} 
}
\details{
The program does not handle ties in the data. We suggest to randomly split the ties before using the program.}
\value{
A table with 3 lines, one for each of the following test: distance from origin, area under the curve (AUC) and restrained adaptive tests. For each test, the value of the statistic and the p-value are given, with the specified number of digits.
}
\references{
	\tabular{l}{\cr
	Chauvel, C, OQuigley, J (2014) Tests for comparing estimated survival functions. \emph{Biometrika} \cr \bold{101}, 3, 535 – 552. \cr
  Chauvel, C (2014). PhD thesis (in French): Processus empiriques pour l'inférence dans le \cr modèle de survie à risques non proportionnels. \cr Université Pierre et Marie Curie - Paris VI.\cr}
}
\author{
Cecile Chauvel
}


\seealso{
\code{\link{standscore}}
\code{\link{plotscore}}
}

\examples{
library(survival)
data(ovarian)

#############################################
# Tests for H_0: beta = 0 for both age and rx covariates

testscore(formula=futime+fustat~age+rx,data=ovarian)

#############################################
# Tests for H_0: beta=  maximum partial likelihood estimator of beta in the Cox model

beta_cox=coxph(Surv(futime,fustat)~ age+rx,data=ovarian)$coeff

testscore(formula=futime+fustat~age+rx,data=ovarian,beta0=beta_cox)

}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
