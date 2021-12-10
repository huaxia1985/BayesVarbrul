library(MCMCpack)
library(mvtnorm)
library(stringr)

pathname <- "~/Desktop/BayesVarbrul" #set working path, e.g. a folder called "BayesVarbrul" on desktop
setwd(pathname)
source("BayesVarbrul.R") #make sure all the files are in the working path

##prepare data
rawdata <- read.csv("rawdata.csv",check.names=F)
tmp <- str_sub(colnames(rawdata),2)
L <- as.numeric(tmp[length(tmp)]) #number of language variables
data <- NULL
factor.name <- c("GENDER","AGE","HH1","HL2","HH2","LL2","LH3","LL3","HH3") #the column names of the social factors in the rawdata that you want to include in the analysis #reference is High exposure Low education G1 (HL1) female
data$X <- t(rawdata[-1,factor.name])
colnames(data$X) <- NULL
data$X <- as.matrix(data$X)
K <- dim(data$X)[2] #number of speakers
data$lingua <- vector("list",L)
for (i in 1:L) {
    data$lingua[[i]] <- rawdata[-1,which(tmp==as.character(i))]
    nl <- sum(tmp==as.character(i))
    if (nl==1) {
    	data$lingua[[i]] <- cbind(data$lingua[[i]],2)  #add all 2 to unrecorded variant, so that it does not affect likelihood
    }
}
data$type <- vector("list",L)
for (i in 1:L) {
	data$type[[i]] <- as.numeric(rawdata[1,which(tmp==as.character(i))])
	nl <- sum(tmp==as.character(i))
    if (nl==1) {
		data$type[[i]] <- c(data$type[[i]],4)
	}
}
data$theta0 <- vector("list",L)
for (i in 1:L) {
	data$theta0[[i]] <- numeric(length(data$type[[i]]))
	tmp2 <- which(data$type[[i]]==1)  #heritage variants has initial frequency = 1 in total
	data$theta0[[i]][tmp2] <- 1/length(tmp2) #if there are more than one heritage variants, then all heritage variants have initial frequency = 1/number of heritage variants
}
n <- max(unlist(data$type))

##start the analysis
BayesVarbrul(data,inter=10^6,interval=10^3,pathname=paste0(pathname,"/result")) #use default for the rest inputs if not sure
#this will take days to finish, depending on the data size
#the results are saved in "result" folder in the working path

##summarize the result
#examine MCMC convergence
post <- read.csv("result/post.csv",sep=" ",header=F)
colnames(post) <- c("sample","interation","logLik","logPrior",factor.name,"G","K","I","U","m")
plot(c(1:dim(post)[1]),post$logPrior+post$logLik,type="l")
burin <- 500
post <- post[-c(1:burin),] # throw away the first few samples as burin
TT <- dim(post)[1]

#plot marginal effect
par(mfrow=c(4,3))

G2 <- (post$HL2+post$HH2+post$LL2)/3
prob <- 1-sum(G2<=0)/TT
den <- density(G2)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Gen2 compared to Gen1",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

G3 <- (post$LH3+post$LL3+post$HH3)/3
prob <- 1-sum(G3<=0)/TT
den <- density(G3)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Gen3 compared to Gen1",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum((G3-G2)<=0)/TT
den <- density(G3-G2)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Gen3 compared to Gen2",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$HH1<=0)/TT
den <- density(post$HH1)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Effect of education in G1",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum((post$HH2-post$HL2)<=0)/TT
den <- density(post$HH2-post$HL2)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <-  1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Effect of education in G2",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum((post$LH3-post$LL3)<=0)/TT
den <- density(post$LH3-post$LL3)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <-  1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Effect of education in G3",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum((post$HL2-post$LL2)<=0)/TT
den <- density(post$HL2-post$LL2)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <-  1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Effect of exposure to Gurindji in G2",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum((post$HH3-post$LH3)<=0)/TT
den <- density(post$HH3-post$LH3)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <-  1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Effect of exposure to Gurindji in G3",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$GENDER<=0)/TT
den <- density(post$GENDER)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <-  1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Male compared to Female",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$AGE<=0)/TT
den <- density(post$AGE)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}

prob <- 1-sum(post$AGE<=0)/TT
den <- density(post$AGE)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="age",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$K<=post$G)/TT
den <- density(post$K/post$G)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <-  1-prob
	tmp <- den$x<=1
}
plot(den,xlab="bias towards Kriol",main=round(prob,2))
abline(v=1)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$I<=post$G)/TT
den <- density(post$I/post$G)
tmp <- den$x>=1
if (prob <= 0.5) {
	prob <-  1-prob
	tmp <- den$x<=1
}
plot(den,xlab="bias towards Innovative",main=round(prob,2))
abline(v=1)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

#calculate Bayes factor
#read in posterior samples of theta and omega
tmp <- read.table("result/omega.csv",header=F,sep=" ",fill=T,stringsAsFactors=F)
idx <- which(is.na(tmp[,K+1]))
omega.post <- vector("list",TT)
for (j in 1:TT) {
	i <- j+burin
	tmp2 <- tmp[(idx[i]+1):(ifelse(i<length(idx),idx[i+1]-1,dim(tmp)[1])),-1]
    omega.post[[j]] <- data.frame(apply(tmp2,2,function (x) as.numeric(x)))
}
theta.post <- vector("list",L)
for (ii in 1:L) {
	filename <- paste0("result/theta",ii,".csv")
	tmp <- read.table(filename,header=F,sep=" ",fill=T,stringsAsFactors=F)
	idx <- which(is.na(tmp[,K+1]))
	theta.post[[ii]] <- vector("list",TT)
	for (j in 1:TT) {
		i <- j+burin
		tmp2 <- tmp[(idx[i]+1):(ifelse(i<length(idx),idx[i+1]-1,dim(tmp)[1])),-1]
    	theta.post[[ii]][[j]] <- data.frame(apply(tmp2,2,function (x) as.numeric(x)))
    }
}

#calculate usage frequency expected by model without effects
nl <- sapply(data$type,function (i) length(i))
theta0.post <- vector("list",L)
N <- 100 #default, fixed to avoid model unidentifiability
for (ii in 1:L) {
	theta0.post[[ii]] <- vector("list",TT-burin+1)
	for (j in 1:TT) {
		if (estimate.r==TRUE) {
		r <- c(post$G[j],post$K[j],post$I[j],post$U[j])
		r_tmp <- r[data$type[[ii]]]
		for (z in 1:n) {
			idx <- which(data$type[[ii]]==z)
			if (length(idx)>1) {
				r_tmp[idx] <- r_tmp[idx]/length(idx)
			}
		}
		m <- post$m[j]
		beta <- as.numeric(post[j,rownames(data$X)])
		a <- N/(N+sum(r_tmp))
		b <- r_tmp/(N+sum(r_tmp))
		tmp <- exp((a-1)*m)*(data$theta0[[ii]]-b/(1-a))+b/(1-a)
		} else {
			tmp <- data$theta0[[ii]]
		}
    	theta0.post[[ii]][[j]] <- matrix(rep(tmp,K),ncol=K,byrow=F)
    }
}

#function to integrate to calculate BF for all factors as total effect
integrad <- function (b,theta0,theta,omega.post,post,X,nl.i,anc,K,TT,factor.name) {
	f <- function (b.i) {
	W <- 0
	for (j in 1:TT) {
		omega <- matrix(as.numeric(unlist(omega.post[[j]])),K,K)
		omega[lower.tri(omega)] <- t(omega)[lower.tri(omega)]
		beta <- as.numeric(post[j,factor.name])
    	theta1 <- theta0[[j]]
		if (nl.i>length(anc)) {
			theta1[-anc,] <- theta1[-anc,]/sum(theta0[[j]][-anc])*matrix(rep(exp(log(sum(theta0[[j]][-anc])/sum(theta0[[j]][anc]))+b.i*beta%*%X),nl.i-length(anc)),nrow=nl.i-length(anc),byrow=T)*sum(theta0[[j]][anc])
			theta1 <- theta1/t(matrix(rep(colSums(theta1),nl.i),K,nl.i))
		}
		prcal0 <- function (theta_hat0,theta_hat,theta_obs,omega,nl.i) {
			V <- matrix(NA,nl.i,nl.i)
			diag(V) <- theta_hat0*(1-theta_hat0)
			for (i in 1:(nl.i-1)) {
				for (j in (i+1):nl.i) {
					V[i,j] <- V[j,i] <- -theta_hat0[i]*theta_hat0[j]
				}
			}
			theta_trans <- theta_hat
			V_trans <- diag(V)
			if (nl.i>2) {
				for (i in 2:(nl.i-1)) {
					tmp <- V[i,(i+1):nl.i]%*%solve(V[(i+1):nl.i,(i+1):nl.i])
					V_trans[i] <- V_trans[i]-tmp%*%V[(i+1):nl.i,i]
					theta_trans[i,] <- theta_hat[i,]+tmp%*%(theta_obs[(i+1):nl.i,]-theta_hat[(i+1):nl.i,])
				}
			}
			lnPr <- 0
			for (i in 2:nl.i) {
				lnPr <- lnPr + dmvnorm(x=theta_obs[i,],mean=theta_trans[i,],sigma=V_trans[i]*omega,log=T)
			}
			lnPr
		}
		lnPr0 <- prcal0(theta0[[j]][,1],theta0[[j]],as.matrix(theta[[j]]),omega,nl.i)
		lnPr1 <- prcal0(theta0[[j]][,1],theta1,as.matrix(theta[[j]]),omega,nl.i)
		tmp <- exp(lnPr1-lnPr0)
		if (!is.finite(tmp) || is.na(tmp)) {tmp <- 0}
		W <- W+tmp
    }
    res <- W/TT
    }
    unlist(sapply(b,function (i) f(i)))
}

#function to integrate to calculate BF for individual factor
integrad.ind <- function (b,theta0,theta,omega.post,post,X,nl.i,anc,K,TT,factor.name,ind.name) {
	f <- function (b.i) {
	W <- 0
	for (j in 1:TT) {
		omega <- matrix(as.numeric(unlist(omega.post[[j]])),K,K)
		omega[lower.tri(omega)] <- t(omega)[lower.tri(omega)]
		beta0 <- beta <- as.numeric(post[j,factor.name])
        #reference in HL1
        Edu1 <- beta[3]  #HH1-HL1
        Edu2 <- beta[5]-beta[4] #HH2-HL2
        Edu3 <- beta[7]-beta[8] #LH3-LL3
        Exp2 <- beta[6]-beta[4] #LL2-HL2
        Exp3 <- beta[9]-beta[7] #HH3-LH3
        G2 <- beta[4] #HL2-HL1
        G3 <- (beta[7]-Edu3-Exp3+beta[8]-Exp3+beta[9]-Edu3)/3 #LH3->HL3 + LL3->HL3 + HH3->HL3
		if (ind.name=="gen") {
			beta[4:6] <- beta[4:6]-G2+b.i*G2
			beta[7:9] <- beta[7:9]-G3+b.i*G3
		}
		if (ind.name=="edu") {
			beta[3] <- b.i*Edu1
			beta[5] <- beta[5]-Edu2+b.i*Edu2
			beta[7] <- beta[7]-Edu3+b.i*Edu3
			beta[9] <- beta[9]-Edu3+b.i*Edu3
		}
		if (ind.name=="exp") {
			beta[6] <- beta[6]-Exp2+b.i*Exp2
			beta[7] <- beta[7]-Exp3+b.i*Exp3
			beta[8] <- beta[8]-Exp3+b.i*Exp3
		}
		if (ind.name=="gender") {
			beta[1] <- b.i*beta[1]
		}
    	theta1 <- theta0[[j]]
		if (nl.i>length(anc)) {
			theta1[-anc,] <- theta1[-anc,]/sum(theta0[[j]][-anc])*matrix(rep(exp(log(sum(theta0[[j]][-anc])/sum(theta0[[j]][anc]))+beta%*%X),nl.i-length(anc)),nrow=nl.i-length(anc),byrow=T)*sum(theta0[[j]][anc])
			theta1 <- theta1/t(matrix(rep(colSums(theta1),nl.i),K,nl.i))
		}
		if (ind.name=="gen") {
			beta0[4:6] <- beta0[4:6]-G2
			beta0[7:9] <- beta0[7:9]-G3
		}
		if (ind.name=="edu") {
			beta0[3] <- 0
			beta0[5] <- beta0[5]-Edu2
			beta0[7] <- beta0[7]-Edu3
			beta0[9] <- beta0[9]-Edu3
		}
		if (ind.name=="exp") {
			beta0[6] <- beta0[6]-Exp2
			beta0[7] <- beta0[7]-Exp3
			beta0[8] <- beta0[8]-Exp3
		}
		if (ind.name=="gender") {
			beta0[1] <- 0
		}
    	theta00 <- theta0[[j]]
		if (nl.i>length(anc)) {
			theta00[-anc,] <- theta00[-anc,]/sum(theta0[[j]][-anc])*matrix(rep(exp(log(sum(theta0[[j]][-anc])/sum(theta0[[j]][anc]))+beta0%*%X),nl.i-length(anc)),nrow=nl.i-length(anc),byrow=T)*sum(theta0[[j]][anc])
			theta00 <- theta00/t(matrix(rep(colSums(theta00),nl.i),K,nl.i))
		}
		prcal0 <- function (theta_hat0,theta_hat,theta_obs,omega,nl.i) {
			V <- matrix(NA,nl.i,nl.i)
			diag(V) <- theta_hat0*(1-theta_hat0)
			for (i in 1:(nl.i-1)) {
				for (j in (i+1):nl.i) {
					V[i,j] <- V[j,i] <- -theta_hat0[i]*theta_hat0[j]
				}
			}
			theta_trans <- theta_hat
			V_trans <- diag(V)
			if (nl.i>2) {
				for (i in 2:(nl.i-1)) {
					tmp <- V[i,(i+1):nl.i]%*%solve(V[(i+1):nl.i,(i+1):nl.i])
					V_trans[i] <- V_trans[i]-tmp%*%V[(i+1):nl.i,i]
					theta_trans[i,] <- theta_hat[i,]+tmp%*%(theta_obs[(i+1):nl.i,]-theta_hat[(i+1):nl.i,])
				}
			}
			lnPr <- 0
			for (i in 2:nl.i) {
				lnPr <- lnPr + dmvnorm(x=theta_obs[i,],mean=theta_trans[i,],sigma=V_trans[i]*omega,log=T)
			}
			lnPr
		}
		lnPr0 <- prcal0(theta0[[j]][,1],theta00,as.matrix(theta[[j]]),omega,nl.i)
		lnPr1 <- prcal0(theta0[[j]][,1],theta1,as.matrix(theta[[j]]),omega,nl.i)
		tmp <- exp(lnPr1-lnPr0)
		if (!is.finite(tmp) || is.na(tmp)) {tmp <- 0}
		W <- W+tmp
    }
    res <- W/TT
    }
    unlist(sapply(b,function (i) f(i)))
}

bf <- numeric(L)
bf.gen <- numeric(L)
bf.edu <- numeric(L)
bf.exp <- numeric(L)
bf.gender <- numeric(L)
for (ii in 1:L) {
	anc <- which(data$type[[ii]]==1)
	bf[ii] <- integrate(integrad,lower=0,upper=100,theta0=theta0.post[[ii]],theta=theta.post[[ii]],omega.post=omega.post,post=post,X=data$X,nl.i=nl[ii],anc=anc,K=K,TT=TT,factor.name=factor.name)$value
	bf.gen[ii] <- integrate(integrad.ind,lower=0,upper=100,theta0=theta0.post[[ii]],theta=theta.post[[ii]],omega.post=omega.post,post=post,X=data$X,nl.i=nl[ii],anc=anc,K=K,TT=TT,factor.name=factor.name,ind.name="gen")$value
	bf.edu[ii] <- integrate(integrad.ind,lower=0,upper=100,theta0=theta0.post[[ii]],theta=theta.post[[ii]],omega.post=omega.post,post=post,X=data$X,nl.i=nl[ii],anc=anc,K=K,TT=TT,factor.name=factor.name,ind.name="edu")$value
	bf.exp[ii] <- integrate(integrad.ind,lower=0,upper=100,theta0=theta0.post[[ii]],theta=theta.post[[ii]],omega.post=omega.post,post=post,X=data$X,nl.i=nl[ii],anc=anc,K=K,TT=TT,factor.name=factor.name,ind.name="exp")$value
	bf.gender[ii] <- integrate(integrad.ind,lower=0,upper=100,theta0=theta0.post[[ii]],theta=theta.post[[ii]],omega.post=omega.post,post=post,X=data$X,nl.i=nl[ii],anc=anc,K=K,TT=TT,factor.name=factor.name,ind.name="gender")$value
}

#plot Bayes factor
var <- read.csv("var.csv")
color <- rep("purple",L)
color[var$type=="PLV"] <- "blue"
color[var$type=="PGN"] <- "orange"
color[var$type=="PGV"] <- "gold"
color[var$type=="CLN"] <- "forestgreen"

plot(log(unlist(bf)),type="h",col=color,xlab="",ylab="log BF",main="Total effect",cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5)
tmp <- which(unlist(bf)>=100)
abline(h=log(100),lty="dashed")
text(x=c(1:L)[tmp],y=log(unlist(bf))[tmp],labels=var$name[tmp],adj=0.5,cex=1)
legend(x=140,y=40,legend=c("Production Lexicon Noun","Production Lexicon Verb","Production Grammar Noun","Production Grammar Verb","Comprehension Lexicon Noun"),fill=c("purple","blue","orange","gold","forestgreen"),cex=1.5)

plot(log(unlist(bf.gen)),type="h",col=color,xlab="",ylab="log BF",main="Generational effect",cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5)
tmp <- which(unlist(bf.gen)>=100)
abline(h=log(100),lty="dashed")
text(x=c(1:L)[tmp],y=log(unlist(bf.gen))[tmp],labels=var$name[tmp],adj=0.5,cex=1)

plot(log(unlist(bf.edu)),type="h",col=color,xlab="",ylab="log BF",main="Effect of education in English",cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5)
tmp <- which(unlist(bf.edu)>=100)
abline(h=log(100),lty="dashed")
text(x=c(1:L)[tmp],y=log(unlist(bf.edu))[tmp],labels=var$name[tmp],adj=0.5,cex=1)

plot(log(unlist(bf.exp)),type="h",col=color,xlab="",ylab="log BF",main="Effect of exposure to Gurindji",cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5)
tmp <- which(unlist(bf.exp)>=100)
abline(h=log(100),lty="dashed")
text(x=c(1:L)[tmp],y=log(unlist(bf.exp))[tmp],labels=var$name[tmp],adj=0.5,cex=1)

plot(log(unlist(bf.gender)),type="h",col=color,xlab="",ylab="log BF",main="Effect of gender",cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5)
tmp <- which(unlist(bf.gender)>=100)
abline(h=log(100),lty="dashed")
text(x=c(1:L)[tmp],y=log(unlist(bf.gender))[tmp],labels=var$name[tmp],adj=0.5,cex=1)

#plot correlation matrix among speakers
omega <- as.matrix(omega.post[[sample.int(TT,1)]])
tmp <- diag(1,K)
diag(tmp) <- diag(omega)^(-1/2)
corr <- tmp%*%omega%*%tmp

a <- heatmap(corr)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colSide <- brewer.pal(max(rawdata[-1,"FAMILY"]),"BrBG")
rowSide <- col_vector[1:max(rawdata[-1,"LOCATION"])]
bk <- c(min(corr),seq(min(corr)+(max(corr[round(corr,2)<1])-min(corr))/30,max(corr[round(corr,2)<1]),by=(max(corr[round(corr,2)<1])-min(corr))/30),1)
color <- c(colorRampPalette(colors = c("blue","white"))(15),colorRampPalette(colors = c("white","orange"))(15),"red")
heatmap(corr,ColSideColors=colSide[rawdata[-1,"FAMILY"]],RowSideColors=rowSide[rawdata[-1,"LOCATION"]],breaks=bk,scale="none",col=color)

#calculate intraclass correlation
intracorr.cal <- function (omega,a) {
	n <- max(a)
	corrvec <- numeric(n)
	for (i in 1:n) {
		tmp <- omega[a==i,a==i]
		if (sum(a==i)>1) {
			corrvec[i] <- mean(tmp[lower.tri(tmp)])
		} else {
			corrvec[i] <- 0
		}
	}
	sum(corrvec)/n/sum(diag(omega))*dim(omega)[1]
}
corr.fam <- numeric(TT)
corr.hou <- numeric(TT)
for (j in 1:TT) {
		omega <- matrix(as.numeric(unlist(omega.post[[j]])),K,K)
		omega[lower.tri(omega)] <- t(omega)[lower.tri(omega)]
		corr.fam[j] <- intracorr.cal(omega,a=rawdata$FAMILY[-1])
		corr.hou[j] <- intracorr.cal(omega,a=rawdata$LOCATION[-1])
	}

#plot intraclass correlation
par(mfrow=c(1,2))
prob <- 1-sum(corr.fam<=0)/TT
den <- density(corr.fam)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Intra-family correlation",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(corr.hou<=0)/TT
den <- density(corr.hou)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Intra-household correlation",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")
