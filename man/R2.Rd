\name{R2}
\alias{R2}
\title{R2 coefficient}
\description{
This function calcultates the \eqn{R^2} coefficient of OQuigley and Flandre (1994) to evaluate the predictive capacity of the proportional hazards model (or Cox model).}
\usage{
R2(formula, data)
}
\arguments{
  \item{formula}{
A formula object or character string with the time and censoring status separated by "+" on the left hand side and the covariates separated by "+" on the right. For instance, if the time name is "Time", the censoring status is "Status" and the covariates are called "Cov1" and "Cov2", the formula is "Time+Status~Cov1+Cov2".
}
  \item{data}{ 
  A data.frame with the data. The censoring status should be 1 for failure and 0 for censoring. No missing data accepted.
  	 }
}

\value{
	\itemize{
\item If one covariate Z is present in the model, the \eqn{R^2} coefficient is 
\deqn{ R^2=1-\frac{\sum(Zi-E_b(Zi))^2}{\sum(Zi-E_0(Zi))^2}, }
where the sums are over the failures. \eqn{E_b(Zi)} is the expectation of \eqn{Z} at the ith failure time under the model of parameter \eqn{b} = the maximum partial likelihood estimator of the regression coefficient. \eqn{E_0(Zi)} is the expectation of \eqn{Z} under the model of parameter 0 at the ith failure time. 

\item If several covariates are present in the model, the \eqn{R^2} coefficient is evaluated as in the previous case except that the covariate Z is replaced by the prognostic index \eqn{b'Z}.
			}
	}
\details{ The program does not handle ties in the data. We suggest to randomly split the ties before using the program.
		}
\references{
OQuigley, J. (2008) \emph{Proportional hazards regression}. Springer New-York. Chapter 12. \cr
OQuigley J, Flandre P. (1994) Predictive capability of proportional hazards regression. \emph{PNAS} \bold{91}, 2310-2314. 
}


\author{Cecile Chauvel}



\examples{
library(survival)
data(ovarian)
R2(futime+fustat~age,data=ovarian)
R2(futime+fustat~age+rx,data=ovarian)
}

%\keyword{ Proportional hazards model }
%\keyword{ Cox Model }
%\keyword{ Survival }
%\keyword{ R2 coefficient }
%\keyword{ Explained Variation }
