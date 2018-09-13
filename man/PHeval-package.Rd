\name{PHeval-package}
\alias{PHeval-package}
\alias{PHeval}
\docType{package}
\title{
Evaluation of the proportional hazards assumption with a standardized score process
}
\description{This package provides tools for the evaluation of the goodness of fit and the predictive capacity of the proportional hazards model.
}
\details{
\tabular{ll}{
Package: \tab PHeval\cr
Type: \tab Package\cr
Version: \tab 0.5.1\cr
Date: \tab 2014-04-20\cr
License: \tab GPL (>=2.0)\cr
}
This package provides functions to evaluate and plot the standardized score process of OQuigley (2003). The plot of this process over the ranks of the failure times gives an indication of the validity of the proportional hazards assumption. A function to evaluate the \eqn{R^2} coefficient of OQuigley and Flandre (1994) is provided. This coefficient is a measure of the predictive ability of the proportional hazards model.
}
\author{
Cecile Chauvel <cecile.chauvel@upmc.fr>
}
\references{
	\tabular{l}{\cr
	
Chauvel, C. and OQuigley, J. (2014) Tests for comparing estimated survival functions. \emph{Biometrika} \cr \bold{101}, 535-552. https://doi.org/10.1093/biomet/asu015 \cr
Chauvel, C. and OQuigley, J. (2017) Survival model construction guided by fit and predictive strength.  \cr \emph{Biometrics} \bold{73}, 483-494. https://doi.org/10.1111/biom.12611 \cr
OQuigley, J. (2008) \emph{Proportional hazards regression}. Springer New-York. \cr		
OQuigley J, Flandre P. (1994) Predictive capability of proportional hazards regression. \emph{PNAS} \bold{91}, 2310-2314. \cr
			}
}

\examples{
library(survival)
S=standscore(futime+fustat~resid.ds+rx,data=ovarian)
plotscore(S,printCB = TRUE)
R2(futime+fustat~resid.ds+rx,data=ovarian)
}
