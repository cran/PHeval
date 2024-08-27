###################################################################
# private functions

E_t=function(t,data,beta0,it,ic)
{ Z=data[,ic][data[,it]>=data[t,it]] 
pi=exp(beta0*Z)/sum(exp(beta0*Z))
sum(Z*pi)
}

EV_t=function(t,data,beta0,it,ic)
{ Z=data[data[,it]>=data[t,it],ic ]
pi=exp(beta0*Z)/sum(exp(beta0*Z))
c(sum(Z*pi),sum(Z^2*pi)-sum(Z*pi)^2)
}

V_t=function(t,data,beta0,it,ic)
{ Z=data[data[,it]>=data[t,it],ic ]
pi=exp(beta0*Z)/sum(exp(beta0*Z))
sum(Z^2*pi)-sum(Z*pi)^2
}


residual=function(t,data,beta0,it,ic)
{data[t,ic]-E_t(t,data,beta0,it,ic)
}

standresidual=function(ic,data,beta0,it,T)
{
EV=sapply(T,FUN=EV_t,data,beta0,it,ic)
k=length(T)-sum(abs(EV[2,])<10^(-15))

cumscore=c(0,cumsum((data[T[1:k],ic]-EV[1,1:k])/sqrt(k*EV[2,1:k])))
names(cumscore)=c(0,data[T[1:k],it])
return(cumscore)
}



EV_t_multiv=function(t,data,beta0,it,ic)
{
Z=as.matrix(data[data[,it]>=data[t,it],ic ])
pi=exp(Z%*%beta0)/sum(exp(Z%*%beta0))
E= t(Z)%*%pi

covpi=function(i1,i2)
{ sum(Z[,i1]*Z[,i2]*pi)-E[i1,1]*E[i2,1]  }
covpi1=function(i1)
{ sapply(1:length(ic),FUN=covpi,i1=i1)  }

V=sapply(1:length(ic),FUN=covpi1)
colnames(V)=colnames(Z)
rownames(V)=colnames(Z)
retu=as.list(c())
retu[[1]]=E
retu[[2]]=V
return(retu)
}



sqrt_inv_mat=function(M)
{X=eigen(M,symmetric=TRUE)
D=diag(X$values^(-1/2))
X$vectors%*%D%*%t(X$vectors)
}

sqrt_mat=function(M)
{X=eigen(M,symmetric=TRUE)
D=diag(X$values^(1/2))
X$vectors%*%D%*%t(X$vectors)
}

confband=function(s,Vbetahat)
{  
Q=sqrt(-log(0.05/(2))/2)
if(dim(s)[2]>1)
	{ k=dim(s)[1]-1
		cb_aux=function(j)
		{ res=matrix(NA,ncol=2,nrow=k+1) 
		 res[,1]=c(0:k)*s[k+1,j]/k-Q*sqrt( sum( (Vbetahat[j,]^2) ) )
		 res[,2]=c(0:k)*s[k+1,j]/k+Q*sqrt( sum( (Vbetahat[j,]^2) ) )
		 return(res)
		}
	res1=lapply(1:dim(s)[2],FUN=cb_aux)		
	}
else
	{k=length(s)-1
	res1=matrix(NA,ncol=2,nrow=k+1) 
	res1[,1]=c(0:k)*s[k+1]/k-Q*rep(sqrt(Vbetahat),k+1)
	res1[,2]=c(0:k)*s[k+1]/k+Q*rep(sqrt(Vbetahat),k+1)
	}
return(res1)
}






standresidual_multiv=function(ic,data,beta0,it,T,globstan)
{
list_i=function(list,i) 
		{list[[i]]
		}
cumsumrow_mat=function(row,mat)
		{cumsum(c(0,mat[row,]))}

EV=sapply(T,FUN=EV_t_multiv,data,beta0,it,ic)
k=length(T)-sum(abs(sapply(EV[2,],det))<10^(-15))
E_unlist=unlist(EV[1,1:k])
EE=matrix(NA,ncol=length(ic),nrow=length(T))
VV=lapply(EV[2,1:k],FUN=sqrt_inv_mat)

for(i in 1 : length(ic)){
	EE[,i]=E_unlist[0:(length(T)-1)*length(ic)+i ]}

inc=function(j)
	{VV[[j]]%*%(as.numeric(data[T[j],ic])-EE[j,])
	}

increment=lapply(1:k,FUN=inc)


cumscore= sapply(1:length(ic),FUN=cumsumrow_mat,sapply(1:k,list_i,list=increment))/sqrt(k)
if(globstan) 	{ SIGMA_min12=sqrt_inv_mat(matrix(rowMeans(sapply(1:k,list_i,list=EV[2,])),ncol=length(ic)))
	cumscore=t(SIGMA_min12%*%t(cumscore))
	}
colnames(cumscore)=names(data)[ic]
rownames(cumscore)=c(0,data[T[1:k],it])
ret=c()
ret[[1]]=cumscore
if(globstan) 	{ ret[[2]]=SIGMA_min12 }
return(ret)
}


integrale=function(x,y){
  idx = 2:length(x)
  return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)}


pvalmax=function(stat,p,n_rep)
{
  Y=rmvnorm(n_rep,mean=c(0,0,0,0),sigma=matrix(c(1,0,sqrt(3)/2,0,0,1,0,sqrt(3)/2,sqrt(3)/2,0,1,0,0,sqrt(3)/2,0,1),ncol=4))
  U1=Y[,1]
  U2=Y[,2]
  A1=Y[,3]
  A2=Y[,4]
  return(1-mean(((U1^2+U2^2)<=stat)*((A1^2+A2^2)<=stat)))
}


###################################################################
# public functions

standscore=function(formula,data,globstan = TRUE, beta0 = 0)
{ 
v=all.vars(as.formula(formula))
f1=as.formula(paste(v[1]," + ",v[2]," ~ ",paste(v[-c(1:2)],collapse=" + "),sep=""))

if(f1!=as.formula(formula)) 
	stop("Incorrect formula type.")
if(!sum(colnames(data)==v[1])) 
	stop(paste("No variable ",v[1]," corresponding to failure times"))
if(!sum(colnames(data)==v[2])) 
	stop(paste("No variable ",v[2]," corresponding to failure indicators"))
check_covariate=function(name)
	{  if(is.na(match(name,colnames(data)))) stop(paste("No covariate",name))
		if(!is.numeric(data[,match(name,colnames(data))])) stop(paste("Covariate ",name," must be numeric"))
	   return(match(name,colnames(data)))
	}
	
ic=sapply(v[-c(1,2)],FUN=check_covariate)
is=which(colnames(data)== v[2])
it=which(colnames(data)== v[1])

if(!is.numeric(data[,it])) stop("The observed times must be numeric values")
if(sum(is.na(data[,c(it,is,ic)]))) stop("Error: presence of missing data")
if(!is.logical(globstan)) 
	stop("Parameter globstan must be logical")
if(length(beta0)==1){if(beta0==0) 
	beta0=rep(0,length(ic)) }
if(!is.numeric(beta0)) 
	stop("Parameter beta0 should be numeric")
if(length(beta0)!=length(ic))
		stop("The length of the parameter beta0 must be the same as the number of covariates.")


if(sum(unique(names(table(data[,is])))==c("0","1"))!=2) 
	{stop("Censoring indicator must be either 0 or 1")}
data=data[order(data[,it]),]
T=which(data[,is]==1)
ret=as.list(c())

if(length(ic)==1) 	
	{res=sapply(ic,FUN=standresidual,data,beta0,it,T)
betahat=coxph(as.formula(paste("Surv(data[,it],data[,is])~",paste(v[-c	(1:2)],collapse=" + "),sep="")),data)$coeff
	Vb=sapply(T,FUN=V_t,data=data,beta0=betahat,it=it,ic=ic)
	Vbetahat=mean(Vb)
	V0=sapply(T,FUN=V_t,data=data,beta0=beta0,it=it,ic=ic)
	Vbeta0=mean(V0)
	cb=confband(res,Vbetahat/Vbeta0)

	ret[[1]]=res
	ret[[2]]=cb
	names(ret)=c("Score","confband")
	}
else 
	{  
	ret=standresidual_multiv(ic,data,beta0,it,T,globstan)
	names(ret)=c("Score")
	if(globstan) {	# no conf. band if no global standardization
					cb=confband(ret[[1]],ret[[2]]) 
					ret=c(ret,cb)
					names(ret)=c("Score","Sigma",paste("confband",v[-c(1,2)], sep=""))
				 }
	}
return(ret)
}




R2=function(formula,data)
{ 
v=all.vars(as.formula(formula))
f1=as.formula(paste(v[1]," + ",v[2]," ~ ",paste(v[-c(1:2)],collapse=" + "),sep=""))

if(f1!=as.formula(formula)) 
	stop("Incorrect formula type.")
if(!sum(colnames(data)==v[1]))
	stop(paste("No variable ",v[1]," corresponding to failure times"))
if(!sum(colnames(data)==v[2])) 
	stop(paste("No variable ",v[2]," corresponding to failure indicators"))
check_covariate=function(name)
	{  if(is.na(match(name,colnames(data)))) stop(paste("No covariate",name))
		if(!is.numeric(data[,match(name,colnames(data))])) stop(paste("Covariate ",name," must be numeric"))
	   return(match(name,colnames(data)))
	}
	
ic=sapply(v[-c(1,2)],FUN=check_covariate)
is=which(colnames(data)== v[2])
it=which(colnames(data)== v[1])

if(!is.numeric(data[,it])) stop("The observed times must be numeric values")
if(sum(is.na(data[,c(it,is,ic)]))) stop("Error: presence of missing data")
if(sum(names(table(data[,is]))==c("0","1"))!=2) 
	{stop("Censoring indicator must be either 0 or 1")}
data=data[order(data[,it]),]
T=which(data[,is]==1)

if(length(ic)==1) 
		 {  b=coxph(Surv(data[,it],data[,is])~data[,ic])$coeff
		    Rb=sapply(T,FUN=residual,data,b,it,ic)^2
		    R0=sapply(T,FUN=residual,data,0,it,ic)^2		
	   	  } 

if(length(ic)>1) 
		{  b=coxph(as.formula(paste("Surv(data[,it],data[,is])~",paste(v[-c	(1:2)],collapse=" + "),sep="")),data)$coeff
	       eta=b%*%t(data[,ic])
           data$eta=t(eta)
		   ie=which(colnames(data)== "eta")
		   Rb=sapply(T,FUN=residual,data,1,it,ie)^2
		   R0=sapply(T,FUN=residual,data,0,it,ie)^2		
	   	  } 
return(1-sum(Rb)/sum(R0))
}


plotscore=function(s,printCB = FALSE,component.num = 1:dim(s[[1]])[2] ,main = "" ,xlab = "Time" ,ylab = "Standardized score" ,ylim)
{
if(dim(s$Score)[2]==1)
	{
	if(missing(ylim)) ylim=range(c(s))
	k=length(s[[1]])-1
	plot(0:k/k,s[[1]],main=main,xlab=xlab,ylab=ylab,ylim=ylim,type="l")
	if((length(s)>1)&&printCB)
		{	
			lines(0:k/k,s[[2]][,1],lty=2)
			lines(0:k/k,s[[2]][,2],lty=2)
		}
	}
else  
	{ if(length(s)==1) printCB=FALSE
	if (max(component.num)>dim(s[[1]])[2] ) stop("Error on component.num.")
		if(missing(ylim)) { ylim=range(c(s[[1]][,component.num]))
		if(printCB) ylim=range(c(ylim,printCB*c(sapply(s,FUN=range)[,2+component.num]))) }
	k=dim(s[[1]])[1]-1
 	matplot(0:k/k,s[[1]][,component.num],main=main,xlab=xlab,ylab=ylab,ylim=ylim,type="l",lty=1)
 		if((length(s)>1)&&printCB)
 			{	
 			plotic=function(j)
 					{
 					lines(0:k/k,s[[j]][,1],lty=2)
					lines(0:k/k,s[[j]][,2],lty=2)
					}
				if(names(s)[2]=="Sigma")
				  	sapply(2+component.num,FUN=plotic)
 			}
	}
	abline(0,0,lty=3)
}



testscore=function(formula,data,beta0=0,n_rep=10^6,digits=5)
{ 
  if(as.integer(n_rep)!=n_rep) stop("Parameter n_rep must be integer")
  if(n_rep<10^5) stop("Parameter n_rep must be higher than 10^5")
  if(as.integer(digits)!=digits) stop("Parameter digits must be integer")
  # tests on other parameters provided in the function standscore  
  
  score=standscore(formula,data,globstan=FALSE,beta0)$Score
  
  test=matrix(NA,ncol=2,nrow=3)
  rownames(test)=c("Distance","AUC","Restrained adaptive")
  colnames(test)=c("Statistic","pvalue")
  p=dim(score)[2]
  if(p==1)
  { 
    score=as.vector(score)
    k=length(score)-1
    J=integrale(0:k/k,score)*sqrt(3)
    p_J=2*(1-pnorm(abs(J)))
    U=score[k+1]
    p_U=2*(1-pnorm(abs(U)))
    T_max=max(abs(J),abs(U))
    p_max=1-pmvnorm(lower=rep(-T_max,2),upper=rep(T_max,2),mean=0,sigma=matrix(c(1,sqrt(3)/2,sqrt(3)/2,1),ncol=2))[1]
  }
  
  else {  
    Jaux=c()
    k=dim(score)[1]-1
    U=sum(score[k+1,]^2)
    for(i in 1 : p)
      Jaux[i]=integrale(0:k/k,score[,i])
    J=3*sum(Jaux^2)
    T_max=max(abs(J),abs(U))
    p_U=1-pchisq(U,df=p)
    p_J=1-pchisq(J,df=p)
    p_max=pvalmax(T_max,p,digits)
  }
  
  test[,1]=c(U,J,T_max)
  test[,2]=c(p_U,p_J,p_max)
  return(round(test,digits))
}