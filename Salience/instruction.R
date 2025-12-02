library(MCMCpack)
library(mvtnorm)
library(stringr)

pathname <- "~/Desktop/Salience" #for example, save the results to a folder called BayesLang_result under Documents folder.
setwd(pathname) #assuming the rawdata is also save in the same folder.

source("BayesVarbrul.R")

##prepare data
rawdata <- read.csv("rawdata.csv",check.names=F)
tmp <- str_sub(colnames(rawdata),2)
var <- as.numeric(tmp)
var <- unique(var[!is.na(var)])
L <- length(var) #number of language variables
data <- NULL
factor.name <- c("GENDER","AGE","HH1","HL2","HH2","LL2","LH3","LL3","HH3") #the column names of the social factors in the rawdata that you want to include in the analysis
ninter <- 2
data$Xind <- t(rawdata[-c(1:(1+ninter)),factor.name])
colnames(data$Xind) <- NULL
data$Xind["AGE",] <- data$Xind["AGE",]-max(data$Xind["AGE",]) #intercept is LH1F aged 100
data$Xind <- as.matrix(data$Xind)
K <- dim(data$Xind)[2] #number of speakers

data$lingua <- vector("list",L)
for (i in 1:L) {
    data$lingua[[i]] <- rawdata[-c(1:(1+ninter)),which(tmp==as.character(var[i]))]
    nl <- sum(tmp==as.character(var[i]))
    if (nl==1) {
    	data$lingua[[i]] <- cbind(data$lingua[[i]],2)  #add all 2 to unrecorded variant, so that it does not affect likelihood
    }
}
data$type <- vector("list",L)
for (i in 1:L) {
	data$type[[i]] <- as.numeric(rawdata[1,which(tmp==as.character(var[i]))])
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

ninter <- 12 #G: Gen1 and Sold, Gen2 and Sold, Gen3 and Sold, Gen1 and Syoung, Gen2 and Syoung, Gen3 and Syoung
			#K: Gen1 and Sold, Gen2 and Sold, Gen3 and Sold, Gen1 and Syoung, Gen2 and Syoung, Gen3 and Syoung
n <- max(unlist(data$type))
data$Xinter <- vector("list",ninter)
Gen2 <- data$Xind["HL2",]+data$Xind["HH2",]+data$Xind["LL2",]
Gen3 <- data$Xind["LH3",]+data$Xind["HH3",]+data$Xind["LL3",]
Gen1 <- numeric(K)+1-Gen2-Gen3
for (j in 1:ninter) {
	data$Xinter[[j]] <- vector("list",L)
	for (i in 1:L) {
		tmp2 <- which(tmp==as.character(var[i]))
		if (j==1) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen1,function (z) log(rawdata[2,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==2),] <- 0
			}
		if (j==2) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen2,function (z) log(rawdata[2,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==2),] <- 0
			}
		if (j==3) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen3,function (z) log(rawdata[2,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==2),] <- 0
			}
		if (j==4) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen1,function (z) log(rawdata[3,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==2),] <- 0
			}
		if (j==5) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen2,function (z) log(rawdata[3,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==2),] <- 0
			}
		if (j==6) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen3,function (z) log(rawdata[3,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==2),] <- 0
			}
		if (j==7) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen1,function (z) log(rawdata[2,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==1),] <- 0
			}
		if (j==8) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen2,function (z) log(rawdata[2,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==1),] <- 0
			}
		if (j==9) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen3,function (z) log(rawdata[2,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==1),] <- 0
			}
		if (j==10) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen1,function (z) log(rawdata[3,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==1),] <- 0
			}
		if (j==11) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen2,function (z) log(rawdata[3,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==1),] <- 0
			}
		if (j==12) {
			data$Xinter[[j]][[i]] <- matrix(as.numeric(sapply(Gen3,function (z) log(rawdata[3,tmp2])*z)),nrow=length(tmp2))
			data$Xinter[[j]][[i]][which(data$type[[i]]==1),] <- 0
			}
	}
}
##start the analysis
BayesVarbrul(data,inter=10^6,interval=10^3,pathname=paste0(pathname,"/result"),estimate.r=F,sd_beta=0.01) #use default for the rest inputs if not sure
#this will take days to finish, depending on the data size

##summarize the result
#examine MCMC convergence
post <- read.csv("~/Desktop/Bayeslangnor/generation_type/L2/post.csv",sep=" ",header=F)
colnames(post) <- c("sample","interation","logLik","logPrior",factor.name,"SoldG1","SoldG2","SoldG3","SyoungG1","SyoungG2","SyoungG3","SoldK1","SoldK2","SoldK3","SyoungK1","SyoungK2","SyoungK3","G","K","m")
plot(c(1:dim(post)[1]),post$logPrior+post$logLik,type="l")
burin <- 10
post <- post[-c(1:burin),] # throw away the first few samples as burin

post2 <- read.csv("~/Desktop/Bayeslangnor/generation_type/L3/post.csv",sep=" ",header=F)
colnames(post2) <- c("sample","interation","logLik","logPrior",factor.name,"SoldG1","SoldG2","SoldG3","SyoungG1","SyoungG2","SyoungG3","SoldK1","SoldK2","SoldK3","SyoungK1","SyoungK2","SyoungK3","G","K","m")
plot(c(1:dim(post2)[1]),post2$logPrior+post2$logLik,type="l")
burin <- 10
post2 <- post2[-c(1:burin),] # throw away the first few samples as burin

post <- rbind(post,post2)

post2 <- read.csv("~/Desktop/Bayeslangnor/generation_type/L1/post.csv",sep=" ",header=F)
colnames(post2) <- c("sample","interation","logLik","logPrior",factor.name,"SoldG1","SoldG2","SoldG3","SyoungG1","SyoungG2","SyoungG3","SoldK1","SoldK2","SoldK3","SyoungK1","SyoungK2","SyoungK3","G","K","m")
plot(c(1:dim(post2)[1]),post2$logPrior+post2$logLik,type="l")
burin <- 10
post2 <- post2[-c(1:burin),] # throw away the first few samples as burin

post <- rbind(post,post2)

TT <- dim(post)[1]
#plot posterior probability for each factor
par(mfrow=c(4,6))
for (i in 1:dim(data$Xind)[1]) {
	plot(c(1:TT),post[,i+4],type="l",xlab=rownames(data$Xind)[i])
}
for (i in 1:length(data$Xinter)) {
	plot(c(1:TT),post[,i+4+dim(data$Xind)[1]],type="l")
}
plot(c(1:TT),post[,i+5+dim(data$Xind)[1]],type="l")
plot(c(1:TT),post[,i+6+dim(data$Xind)[1]],type="l")
plot(c(1:TT),post[,i+7+dim(data$Xind)[1]],type="l")

#effetive size
library(coda)
effectiveSize(as.mcmc(post[,-c(1:4)]))

#test chain convergence
chains <- list(as.mcmc(post[1:990,-c(1:4,26:28)]),as.mcmc(post[c(1:990)+990,-c(1:4,26:28)]),as.mcmc(post[c(1:990)+990*2,-c(1:4,26:28)]))
gelman.diag(chains)
geweke.diag(chains)
heidel.diag(post[,-c(1:4,26:28)])

#plot marginal effect
par(mfrow=c(4,6))

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
plot(den,xlab="age",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SoldG1<=0)/TT
den <- density(post$SoldG1)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with elder in G variants in G1",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SoldG2<=0)/TT
den <- density(post$SoldG2)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with elder in G variants in G2",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SoldG3<=0)/TT
den <- density(post$SoldG3)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with elder in G variants in G3",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SoldK1<=0)/TT
den <- density(post$SoldK1)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with elder in K variants in G1",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SoldK2<=0)/TT
den <- density(post$SoldK2)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with elder in K variants in G2",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SoldK3<=0)/TT
den <- density(post$SoldK3)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with elder in K variants in G3",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SyoungG1<=0)/TT
den <- density(post$SyoungG1)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with young people in G variants in G1",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SyoungG2<=0)/TT
den <- density(post$SyoungG2)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with young people in G variants in G2",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SyoungG3<=0)/TT
den <- density(post$SyoungG3)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with young people in G variants in G3",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SyoungK1<=0)/TT
den <- density(post$SyoungK1)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with young people in K variants in Gen1",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SyoungK2<=0)/TT
den <- density(post$SyoungK2)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with young people in K variants in Gen2",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")

prob <- 1-sum(post$SyoungK3<=0)/TT
den <- density(post$SyoungK3)
tmp <- den$x>=0
if (prob <= 0.5) {
	prob <- 1-prob
	tmp <- den$x<=0
}
plot(den,xlab="Saliency association with young people in K variants in Gen3",main=round(prob,2))
abline(v=0)
polygon(c(den$x[tmp],rev(den$x[tmp])),c(den$y[tmp],numeric(sum(tmp))),col="lightgrey")
