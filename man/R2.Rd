\name{R2}
\alias{R2}
\title{R2 coefficient}
\description{
This function calcultates the R2 coefficient of OQuigley and Flandre (1994) to evaluate the predictive capacity of the proportional hazards model (or Cox model).}
\usage{
R2(formula, data)
}
\arguments{
  \item{formula}{
A formula object or character string with the time and censoring status separated by "+" on the left hand side and the covariates separated by "+" on the right. For instance, if the time name is "Time", the censoring status is "Status" and the covariates are called "Cov1" and "Cov2", the formula is "Time+Status~Cov1+Cov2".
}
  \item{data}{ 
  A data.frame with the data. The censoring status should by 1 for failure and 0 for censoring. No missing data accepted.
  	 }
}

\value{
	\itemize{
\item If one covariate Z is present in the model, the R2 coefficient is \cr
R2=1-sum[(Zi-Eb(Zi))^2]/sum[(Zi-E0(Zi))^2], \cr
where the sums are over the failures, Eb(Zi) is the expectation of Z under the model of parameter the maximum partial likelihood estimator at the ith failure time and E0(Zi) is the expectation of Z under the model of parameter 0 at the ith failure time. 

\item If several covariates are present in the model, the R2 coefficient is \cr
R2=1-sum[(PIi-Eb(PIi))^2]/sum[(PIi-E0(PIi))^2], \cr
where PI=b'Z is the prognostic index evaluated with b = the maximum partial likelihood estimator and Z is the vector of covariates, the sums are over the failures, Eb(PIi) is the expectation of PI under the model of parameter b at the ith failure time and E0(PIi) is the expectation of Z under the model of parameter 0 at the ith failure time.
			}
	}
\details{ The program does not handle ties in the data. We suggest to randomly split the ties before using the program.
		}
\references{
OQuigley J, Flandre P. (1994) Predictive capability of proportional hazards regression. PNAS 91, 2310-2314. \cr
OQuigley, J. (2008) Proportional hazards regression. Springer New-York. Chapter 12.}
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
