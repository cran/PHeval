\name{plotscore}
\alias{plotscore}
\title{ Plot of the standardized score process}
\description{
The function plots the standardized score process. This process is a list resulting of a call to the function standscore. }
\usage{
plotscore(s, printCB = FALSE , component.num = 1:dim(s[[1]])[2], 
main = "" , xlab = "Time", ylab = "Standardized score", ylim)
}

\arguments{
  \item{s}{ A list resulting from a call to the function standscore which corresponds to the standardized score process and, possibly, its confidence bands. }
  \item{printCB}{Set printCB = TRUE for plotting the confidence bands. Default is printCB = FALSE.
}
  \item{component.num}{ A vector of length equals to or lower than the number of covariates. It indicates which components of the process to plot. For example, set component.num=c(2,3) for plotting the second and third components. By default, all components are plotted.}

  \item{main,xlab,ylab,ylim}{ Usual arguments for plotting. See help(plot). }
}

\author{
Cecile Chauvel}

\seealso{
\code{\link{standscore}}
%\code{\link{{par}}
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
score3=standscore(futime+fustat~age+rx,data=ovarian,globstan=FALSE)
plotscore(score3)	

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
