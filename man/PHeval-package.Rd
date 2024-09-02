\name{PHeval-package}
\alias{PHeval-package}
\alias{PHeval}
\docType{package}
\title{
Evaluation of the proportional hazards assumption with a standardized score process
}
\description{This package provides tools for the evaluation of the goodness of fit of the proportional hazards model.
}
\details{
\tabular{ll}{
Package: \tab PHeval\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2024-09-02\cr
License: \tab GPL (>=2.0)\cr
}
This package provides functions to compute and plot the standardized score process of OQuigley (2008). The mathematical properties of this process are studied in Chauvel (2014). The restrained adaptive test introduced in Chauvel and OQuigley (2014) is implemented. The plot of this process over the ranks of the failure times gives an indication of the validity of the proportional hazards assumption. A function to evaluate the \eqn{R^2} coefficient of OQuigley and Flandre (1994) is provided. This coefficient is a measure of the predictive ability of the proportional hazards model.
}
\author{
Cecile Chauvel <cecile.chauvel@univ-lyon1.fr>
}
\references{
	\tabular{l}{\cr
	 Chauvel, C (2014). PhD thesis (in French): Processus empiriques pour l'inférence dans le \cr modèle de survie à risques non proportionnels. \cr Université Pierre et Marie Curie - Paris VI.\cr
	Chauvel, C, OQuigley, J (2014). Tests for comparing estimated survival functions. \cr \emph{Biometrika} 
  \bold{101}, 3, 535 – 552. \cr
OQuigley, J (2008). \emph{Proportional hazards regression}. Springer New-York. \cr
OQuigley, J, Flandre, P (1994). Predictive capability of proportional hazards regression. \emph{PNAS} \bold{91}, 2310-2314. \cr
			}
}

\examples{
library(survival)
S=standscore(futime+fustat~resid.ds+rx,data=ovarian)
plotscore(S,printCB = TRUE)
testscore(futime+fustat~resid.ds+rx,data=ovarian)
R2(futime+fustat~resid.ds+rx,data=ovarian)
}
