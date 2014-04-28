\name{PHeval-package}
\alias{PHeval-package}
\alias{PHeval}
\docType{package}
\title{
Evaluation of the proportional hazards assumption with a standardized score process
}
\description{This package enables to evaluate the standardized score process of OQuigley (2003). A function to evaluate the predictive ability of the proportional hazards model with the R2 coefficient of OQuigley and Flandre (1994) is provided. 
}
\details{
\tabular{ll}{
Package: \tab PHeval\cr
Type: \tab Package\cr
Version: \tab 0.5\cr
Date: \tab 2014-04-20\cr
License: \tab GPL (>=2.0)\cr
}
This package enables to evaluate the standardized score process of OQuigley (2003). The plot of this process over the ranks of the failure times gives an indication of the validity of the proportional hazards assumption. A function to evaluate the predictive ability of the proportional hazards model with the R2 coefficient of OQuigley and Flandre (1994) is provided. 
}
\author{
Cecile Chauvel <cecile.chauvel@upmc.fr>
}
\references{
	\tabular{l}{
OQuigley, J. (2008) Proportional hazards regression. Springer New-York. \cr		
OQuigley, J. (2003) Khmaladze-Type Graphical Evaluation of the Proportional Hazards Assumption.  \cr Biometrika  90, 577-594 \cr 
OQuigley J, Flandre P. (1994) Predictive capability of proportional hazards regression. PNAS 91. \cr
			}
}

\examples{
library(survival)
S=standscore(futime+fustat~resid.ds+rx,data=ovarian)
plotscore(S,printCB = TRUE)
R2(futime+fustat~resid.ds+rx,data=ovarian)
}
