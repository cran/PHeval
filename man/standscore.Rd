\name{standscore}
\alias{standscore}
\title{Standardized Score Process}
\description{
This function evaluates the standardized score process. The process helps evaluating the goodness of fit of the proportional hazards model and visualizing the shape of time-dependent effects.}
\usage{
standscore(formula, data, globstan = TRUE, beta0 = 0)
}
\arguments{
  \item{formula}{
A formula object or character string with the time and censoring status separated by "+" on the left hand side and the covariates separated by "+" on the right. For instance, if the time name is "Time", the censoring status is "Status" and the covariates are "Cov1" and "Cov2", the formula is "Time+Status~Cov1+Cov2".
}
  \item{data}{ 
  A data.frame with the data. The censoring status should be 1 for failure and 0 for censoring. No missing data accepted. } 
 \item{globstan}{
   	With one covariate in the model, globstan has no effect. With several covariates, 
   	\itemize{ 
   		\item if globstan = TRUE (default) a global standardization by the matrix Sigma is applied to the process. Sigma is the estimator of the variance-covariance matrix between the covariates to the power of - 1 / 2. With this sandardization, each component of the process represents the cumulative regression effect of each covariates.
   		\item if globstan = FALSE, no standardization is applied. In this case, the components of the process are dependent and do not reflect the shapes of the cumulative effects. The confidence bands are not given. The use of globstan = FALSE is aimed at performing tests of the value of the regression coefficients (null hypothesis : regression parameter = beta0).
   	}}
  \item{beta0}{
A vector of parameters to evaluate the process (by default, parameters set to 0). Its length is the number of covariates. Each value corresponds to the regression coefficient for a covariate, in the same order as appearing in formula. 
}  

}
\details{
The program does not handle ties in the data. We suggest to randomly split the ties before using the program.}
\value{
  \item{Score}{A vector or matrix with the value of the standardized score process. Each row corresponds to a failure time, each column to a covariate.
  	}
  \item{Sigma }{The matrix used for the standardization of the process. Sigma is the estimator of the variance-covariance matrix between the covariates to the power of - 1 / 2. This value is present only with multiple covariates and globstan = TRUE.
  	}
  \item{confbandCOV}{A matrix with the confidence bands of the process for a constant regression effect associated with the covariate named COV. Each row corresponds to a failure time. The first column is the lower band and the second column is the upper band. This value is present with one covariate or with multiple covariates and globstan = TRUE.  }	
}

\references{
	\tabular{l}{\cr
	Chauvel, C (2014). PhD thesis (in French): Processus empiriques pour l'inf\'erence dans le mod\`ele de survie \`a risques non proportionnels. Universit\'e Pierre et Marie Curie - Paris VI.\cr
	Chauvel, C, OQuigley, J (2014) Tests for comparing estimated survival functions. \cr \emph{Biometrika} 
} 
}

\author{
Cecile Chauvel
}


\seealso{
\code{\link{plotscore}}
}
\examples{
library(survival)
data(ovarian)

#############################################
# Evaluation and plot of the standardized score process at parameter beta0 = 0

score1=standscore(futime+fustat~age+rx,data=ovarian)
plotscore(score1,printCB=TRUE)


#############################################
# Evaluation of the standardized score process at parameter
# beta0 = maximum partial likelihood estimator of beta

beta=coxph(Surv(futime,fustat)~age+rx,data=ovarian)$coeff
score2=standscore(futime+fustat~age+rx,data=ovarian,beta0=beta)

# Separated plots for each regression effect

par(mfrow=c(1,2))
plotscore(score2,printCB=TRUE,component.num=1,main="age")
plotscore(score2,printCB=TRUE,component.num=2,main="rx")

#############################################
# Evaluation and plot of the standardized score process at parameter beta0 = 0
# without global standardization 

fo="futime+fustat~age+rx"
score3=standscore(fo,data=ovarian,globstan=FALSE)
plotscore(score3)
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
