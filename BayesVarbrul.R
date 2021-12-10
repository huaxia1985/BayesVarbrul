##Parameters to estimate: 
#m = number of memory updates. Its prior is a gamma distribution with mean=mu_m, so alpha=beta*mu_m, beta=rate_m.
#r = a vector with each element models the bias to use a variant type. If estimate.r=F, bias is not modelled. r can be thought of as the parameters for a Dirichlet distribution, so to remove their interdependence in the proposal, r is transformed into sr (the sum of r) and rr=r/sr except for the last element in vector r, which is rU for unrecorded types. Then rU=sr(1-sum(rr)). Prior for each r is a gamma distribution with mean=1, so alpha=beta=rate_r.
#theta = a list, with each cell corresponding to a variable. Each cell contains a matrix of usage frequency of each variant by each speaker. Rows are variants and columns are speakers.
#omega = covariance matrix among speakers.
#beta = a vector of coefficients of each social factor, descrbing the amount of increase in the relative frequency of using non-heritage variants to heritage variants when the social factor (if continous) increases, or when a speaker belong to a level of the social factor compared to the reference level. The prior for each beta is normal with mean=0, sd=sd_beta. 

##Required Inputs:
#data
#data$X a matrix records social factors (row) of each speaker (col)
#data$lingua is a list, with each cell corresponding to a variable. Each cell contains a matrix, in which rows are speakers and columns are variants. Binary data is recorded as 0 and 1 for not use and use. Frequency data is recorded as integers, the number of times a variant is used.
#data$type is a list, with each cell corresponding to a variable. Each cell contains a vector, in which each element is the type of each variant in the same order as the variants in data$lingua. Heritage type = 1; Unrecorded type = 4. Unrecorded type is added for each variable with only one variant, such as comprehension variable.
#data$theta0 is a list, with each cell corresponding to a variable. Each cell contains a vector, in whcih each element is the initial usage frequency of each variant in the same order as the variants in data$lingua.
#inter = number of MCMC interations
#interval = number of intervals to write MCMCs to a file
#pathname = the path of the folder to save results
#mu_m = mean of gamma prior for m, also set as the starting value of m
#estimate.r = if bias is modelled
#eta_sr = #tuning parameter for sr
#eta_lr = #tuning parameter for each lr
#eta_m = #tuning parameter for m
#eta_beta = #tuning parameter for each beta
#rate_r = #beta of gamma prior for each r, this give a very large variance to the prior
#rate_m = #beta of gamma prior for m, this give a very large variance to the prior
#sd_beta = #s.d. of normal prior for each beta, this give a very large variance to the prior

##Outputs are in saved files:
#post.csv = each row is a sampled MCMC after each interval, recording interation #,log-likelihood,log-prior,beta,r,m
#thetai.csv = the MCMC sample of theta for the ith language variable. Each sample is a matrix, with rows are variants and columns are speakers.
#omega.csv = the MCMC sample of the covariance matrix among speakers. Each sample is a KbyK matrix, where K is the number of speakers.

BayesVarbrul <- function (data,inter,interval,pathname,mu_m=100,estimate.r=TRUE,eta_sr=0.01,eta_lr=0.1,eta_m=0.1,eta_beta=0.01,rate_r=0.01,rate_m=0.01,sd_beta=100) {
    #check if file exists
    if (file.exists(paste0(pathname,"/post.csv"))) {
        stop("post.csv file exists already, please double check")
    }
    if (file.exists(paste0(pathname,"/omega.csv"))) {
        stop("omega.csv file exists already, please double check")
    }
    L <- length(data$type) #L number of loci
    for (i in 1:L) {
        filename <- paste0(pathname,"/theta",i,".csv")
        if (file.exists(filename)) {
            stop(paste0(filename, "exists already, please double check"))
        }
    }
	#initializing
	K <- dim(data$X)[2] #K number of speakers
	P <- dim(data$X)[1] #P number of predictors
	n <- max(unlist(data$type)) #n number of variant types
	nl <- sapply(data$type,function (i) length(i)) #nl number of variants for each VARIABLE
	rho <- K+sum(nl)
	sr <- 1 #starting value of the sum of r
	rr <- rep(1/n,n-1) #starting value of r/sr
	r <- c(rr,1-sum(rr))*sr
	N <- 100 #memory size, default as a fixed value for model identifiability
	m <- mu_m 
	R <- K*diag(1,K)
	lnL_list <- numeric(L) #list for log-likelihood, with each element corresponding to each language variable
	lnPr_list <- numeric(L) #llist for log-prior, with each element corresponding to each language variable
	theta_hat0 <- vector("list",L) #usage frequency expected from null model
	theta_hat <- vector("list",L) #usage frequency expected from null model and effects of social factors
	theta <- vector("list",L) 
	omega <- riwish(rho,R)
	S <- matrix(0,K,K)
	beta <- numeric(P)
	for (i in 1:L) {
		if (estimate.r==TRUE) {
		r_tmp <- r[data$type[[i]]]
		for (z in 1:n) {
			idx <- which(data$type[[i]]==z)
			if (length(idx)>1) {
				r_tmp[idx] <- r_tmp[idx]/length(idx)
			}
		} 
		a <- N/(N+sum(r_tmp))
		b <- r_tmp/(N+sum(r_tmp))
		theta_hat0[[i]] <- exp((a-1)*m)*(data$theta0[[i]]-b/(1-a))+b/(1-a)
		} else {
			theta_hat0[[i]] <- data$theta0[[i]]
		}
		theta_hat[[i]] <- matrix(rep(theta_hat0[[i]],K),ncol=K,byrow=F)
		anc <- which(data$type[[i]]==1) #anc tells which variants are heritage type.
		if (nl[i]>length(anc)) {
			theta_hat[[i]][-anc,] <- theta_hat[[i]][-anc,]/sum(theta_hat0[[i]][-anc])*matrix(rep(exp(log(sum(theta_hat0[[i]][-anc])/sum(theta_hat0[[i]][anc]))+beta%*%data$X),nl[i]-length(anc)),nrow=nl[i]-length(anc),byrow=T)*sum(theta_hat0[[i]][anc])
			theta_hat[[i]] <- theta_hat[[i]]/t(matrix(rep(colSums(theta_hat[[i]]),nl[i]),K,nl[i]))
		}
		theta[[i]] <- theta_hat[[i]]
		lnL_list[i] <- likcal(data$lingua[[i]],theta[[i]])
		lnPr_list[i] <- prcal(theta_hat0[[i]],theta_hat[[i]],theta[[i]],omega)
	}
	lnL <- sum(lnL_list) #starting value of log likelihood
	lnPr <- sum(lnPr_list)+sum(dgamma(r,shape=rate_r,rate=rate_r,log=T))+sum(dnorm(beta,0,sd_beta,log=T))+dgamma(m,shape=rate_m,rate=rate_m,log=T) #starting value of log prior
	#start MCMC
	for (tt in 1:inter) {
		C <- chol(omega)
		#updating theta
		for (i in 1:L) {
			theta_new <- theta[[i]]+t(t(C)%*%matrix(rnorm(K*nl[i],mean=0,sd=1),K,nl[i]))
			theta_new[theta_new<0] <- 0
			theta_new[theta_new>1] <- 1
			theta_new <- theta_new/t(matrix(rep(colSums(theta_new),nl[i]),K,nl[i]))
			while (is.na(sum(theta_new))) {
				theta_new <- theta[[i]]+t(t(C)%*%matrix(rnorm(K*nl[i],mean=0,sd=0.01),K,nl[i]))
				theta_new[theta_new<0] <- 0
				theta_new[theta_new>1] <- 1
				theta_new <- theta_new/t(matrix(rep(colSums(theta_new),nl[i]),K,nl[i]))
			}
			lnL_tmp <- likcal(data$lingua[[i]],theta_new)
			lnL_new <- lnL-lnL_list[i]+lnL_tmp 
			lnPr_tmp <- prcal(theta_hat0[[i]],theta_hat[[i]],theta_new,omega)
			lnPr_new <- lnPr-lnPr_list[i]+lnPr_tmp
			accept <- min(1,exp(lnL_new+lnPr_new-lnL-lnPr))
			if (runif(1)<accept) {
				lnL <- lnL_new
				lnPr <- lnPr_new
				lnL_list[i] <- lnL_tmp
				lnPr_list[i] <- lnPr_tmp
				theta[[i]] <- theta_new
			}
		}
		#updating beta
		for (j in 1:P) {
			beta_new <- beta
			u <- runif(1,min=-eta_beta,max=eta_beta)
			beta_new[j] <- beta[j]+u
			lnPr_list_new <- numeric(L)
			S_new <- matrix(0,K,K)
			theta_hat_new <- vector("list",L)
			for (i in 1:L) {
				theta_hat_new[[i]] <- matrix(rep(theta_hat0[[i]],K),ncol=K,byrow=F)
				anc <- which(data$type[[i]]==1)
				if (nl[i]>length(anc)) {
					theta_hat_new[[i]][-anc,] <- theta_hat_new[[i]][-anc,]/sum(theta_hat0[[i]][-anc])*matrix(rep(exp(log(sum(theta_hat0[[i]][-anc])/sum(theta_hat0[[i]][anc]))+beta_new%*%data$X),nl[i]-length(anc)),nrow=nl[i]-length(anc),byrow=T)*sum(theta_hat0[[i]][anc])
					theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
				}
				lnPr_list_new[i] <- prcal(theta_hat0[[i]],theta_hat_new[[i]],theta[[i]],omega)
				S_new <- S_new+scal(theta_hat0[[i]],theta_hat_new[[i]],theta[[i]])
			}
			lnPr_new <- lnPr-sum(lnPr_list)+sum(lnPr_list_new)-dnorm(beta[j],0,sd_beta,log=T)+dnorm(beta_new[j],0,sd_beta,log=T)
			accept <- min(1,exp(lnPr_new-lnPr))
			if (runif(1)<accept) {
				lnPr <- lnPr_new
				lnPr_list <- lnPr_list_new
				theta_hat[[i]] <- theta_hat_new[[i]]
				beta <- beta_new
				S <- S_new
			}
		}
		#updating r
		if (estimate.r==TRUE) {
		#updating sr
		u <- runif(1)
		sr_new <- sr*exp(eta_sr*(u-0.5))
		r_new <- c(rr,1-sum(rr))*sr_new
		lnPr_list_new <- numeric(L)
		S_new <- matrix(0,K,K)
		theta_hat_new <- vector("list",L)
		theta_hat0_new <- vector("list",L)
		for (i in 1:L) {
			r_tmp <- r_new[data$type[[i]]]
			for (z in 1:n) {
			idx <- which(data$type[[i]]==z)
			if (length(idx)>1) {
				r_tmp[idx] <- r_tmp[idx]/length(idx)
			}
			} 
			a <- N/(N+sum(r_tmp))
			b <- r_tmp/(N+sum(r_tmp))
			theta_hat0_new[[i]] <- exp((a-1)*m)*(data$theta0[[i]]-b/(1-a))+b/(1-a) 
			theta_hat_new[[i]] <- matrix(rep(theta_hat0_new[[i]],K),ncol=K,byrow=F)
			anc <- which(data$type[[i]]==1)
			if (nl[i]>length(anc)) {
				theta_hat_new[[i]][-anc,] <- theta_hat_new[[i]][-anc,]/sum(theta_hat0_new[[i]][-anc])*matrix(rep(exp(log(sum(theta_hat0_new[[i]][-anc])/sum(theta_hat0_new[[i]][anc]))+beta%*%data$X),nl[i]-length(anc)),nrow=nl[i]-length(anc),byrow=T)*sum(theta_hat0_new[[i]][anc])
				theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
			}
			lnPr_list_new[i] <- prcal(theta_hat0_new[[i]],theta_hat_new[[i]],theta[[i]],omega)
			S_new <- S_new+scal(theta_hat0_new[[i]],theta_hat_new[[i]],theta[[i]])
		}
		lnPr_new <- lnPr-sum(lnPr_list)+sum(lnPr_list_new)-sum(dgamma(r,shape=rate_r,rate=rate_r,log=T))+sum(dgamma(r_new,shape=rate_r,rate=rate_r,log=T))
		accept <- min(1,exp(lnPr_new-lnPr)*exp(eta_sr*(u-0.5)))
		if (runif(1)<accept) {
			lnPr <- lnPr_new
			lnPr_list <- lnPr_list_new
			theta_hat[[i]] <- theta_hat_new[[i]]
			theta_hat0[[i]] <- theta_hat0_new[[i]]
			sr <- sr_new
			r <- r_new
			S <- S_new
		}
		#updating lr
		for (j in 1:(n-1)) {
			lrr <- log(rr[j]/(1-rr[j]))
			u <- runif(1,min=-eta_lr,max=eta_lr)
			lrr_new <- lrr+u
			rr_new <- rr
			rr_new[j] <- 1/(1+exp(-lrr_new))
			r_new <- c(rr_new,1-sum(rr_new))*sr
			lnPr_list_new <- numeric(L)
			S_new <- matrix(0,K,K)
			theta_hat_new <- vector("list",L)
			theta_hat0_new <- vector("list",L)
			for (i in 1:L) {
				r_tmp <- r_new[data$type[[i]]]
				for (z in 1:n) {
				idx <- which(data$type[[i]]==z)
				if (length(idx)>1) {
					r_tmp[idx] <- r_tmp[idx]/length(idx)
				}
				} 
				a <- N/(N+sum(r_tmp))
				b <- r_tmp/(N+sum(r_tmp))
				theta_hat0_new[[i]] <- exp((a-1)*m)*(data$theta0[[i]]-b/(1-a))+b/(1-a) 
				theta_hat_new[[i]] <- matrix(rep(theta_hat0_new[[i]],K),ncol=K,byrow=F)
				anc <- which(data$type[[i]]==1)
				if (nl[i]>length(anc)) {
					theta_hat_new[[i]][-anc,] <- theta_hat_new[[i]][-anc,]/sum(theta_hat0_new[[i]][-anc])*matrix(rep(exp(log(sum(theta_hat0_new[[i]][-anc])/sum(theta_hat0_new[[i]][anc]))+beta%*%data$X),nl[i]-length(anc)),nrow=nl[i]-length(anc),byrow=T)*sum(theta_hat0_new[[i]][anc])
					theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
				}
				lnPr_list_new[i] <- prcal(theta_hat0_new[[i]],theta_hat_new[[i]],theta[[i]],omega)
				S_new <- S_new+scal(theta_hat0_new[[i]],theta_hat_new[[i]],theta[[i]])
			}
			lnPr_new <- lnPr-sum(lnPr_list)+sum(lnPr_list_new)-dgamma(r[j],shape=rate_r,rate=rate_r,log=T)+dgamma(r_new[j],shape=rate_r,rate=rate_r,log=T)
			#jacobian adjustment is f(lr) = 1/(1+exp(-lr)) * (1-1/(1+exp(-lr))) 
			accept <- min(1,exp(lnPr_new-lnPr) * (1/(1+exp(-rr[j]))*(1-1/(1+exp(-rr[j])))) / (1/(1+exp(-rr_new[j]))*(1-1/(1+exp(-rr_new[j])))))
			if (runif(1)<accept) {
				lnPr <- lnPr_new
				lnPr_list <- lnPr_list_new
				theta_hat[[i]] <- theta_hat_new[[i]]
				theta_hat0[[i]] <- theta_hat0_new[[i]]
				rr <- rr_new
				r <- r_new
				S <- S_new
			}
		}
		}
		#updating m
		u <- runif(1)
		m_new <- m*exp(eta_m*(u-0.5))
		lnPr_list_new <- numeric(L)
		S_new <- matrix(0,K,K)
		theta_hat_new <- vector("list",L)
		theta_hat0_new <- vector("list",L)
		for (i in 1:L) {
			if (estimate.r==TRUE) {
			r_tmp <- r[data$type[[i]]]
			for (z in 1:n) {
				idx <- which(data$type[[i]]==z)
				if (length(idx)>1) {
					r_tmp[idx] <- r_tmp[idx]/length(idx)
				}
			} 
			a <- N/(N+sum(r_tmp))
			b <- r_tmp/(N+sum(r_tmp))
			theta_hat0_new[[i]] <- exp((a-1)*m_new)*(data$theta0[[i]]-b/(1-a))+b/(1-a) 
			} else {
				theta_hat0_new[[i]] <- data$theta0[[i]]
			}
			theta_hat_new[[i]] <- matrix(rep(theta_hat0_new[[i]],K),ncol=K,byrow=F)
			anc <- which(data$type[[i]]==1)
			if (nl[i]>length(anc)) {
				theta_hat_new[[i]][-anc,] <- theta_hat_new[[i]][-anc,]/sum(theta_hat0_new[[i]][-anc])*matrix(rep(exp(log(sum(theta_hat0_new[[i]][-anc])/sum(theta_hat0_new[[i]][anc]))+beta%*%data$X),nl[i]-length(anc)),nrow=nl[i]-length(anc),byrow=T)*sum(theta_hat0_new[[i]][anc])
				theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
			}
			lnPr_list_new[i] <- prcal(theta_hat0_new[[i]],theta_hat_new[[i]],theta[[i]],omega)
			S_new <- S_new+scal(theta_hat0_new[[i]],theta_hat_new[[i]],theta[[i]])
		}
		lnPr_new <- lnPr-sum(lnPr_list)+sum(lnPr_list_new)-dgamma(m,shape=rate_m,rate=rate_m,log=T)+dgamma(m_new,shape=rate_m,rate=rate_m,log=T)
		accept <- min(1,exp(lnPr_new-lnPr)*exp(eta_m*(u-0.5)))
		if (runif(1)<accept) {
			lnPr <- lnPr_new
			lnPr_list <- lnPr_list_new
			theta_hat[[i]] <- theta_hat_new[[i]]
			theta_hat0[[i]] <- theta_hat0_new[[i]]
			m <- m_new
			S <- S_new
		}
		#updating omega
		omega <- riwish(rho,R+S)
		for (i in 1:L) {
			lnPr_list[i] <- prcal(theta_hat0[[i]],theta_hat[[i]],theta[[i]],omega)
		}
		lnPr <- sum(lnPr_list)+sum(dgamma(r,shape=rate_r,rate=rate_r,log=T))+sum(dnorm(beta,0,sd_beta,log=T))+dgamma(m,shape=rate_m,rate=rate_m,log=T)
			
	if (tt%%interval==0) {
       print(paste0(tt,"MCMC iteractions have been done"))
       write.table(t(c(tt,lnL,lnPr,beta,r,m)),file=paste0(pathname,"/post.csv"),append=T,row.names=T,col.names=F)
		for (i in 1:L) {
			filename <- paste0(pathname,"/theta",i,".csv") 
			write.table(theta[[i]],file=filename,append=T)
		}
		write.table(omega,file=paste0(pathname,"/omega.csv"),append=T)
	}
	}
}

##Internal functions
scal <- function (theta_hat0,theta_hat,theta) {
	nl <- length(theta_hat0)
	K <- dim(theta)[2]
	V <- matrix(NA,nl,nl)
	diag(V) <- theta_hat0*(1-theta_hat0)
	for (i in 1:(nl-1)) {
		for (j in (i+1):nl) {
			V[i,j] <- V[j,i] <- -theta_hat0[i]*theta_hat0[j]
		}
	}
	theta_trans <- theta_hat
	V_trans <- diag(V)
	S <- 0
	if (nl>2) {
	for (i in 2:(nl-1)) {
		tmp <- V[i,(i+1):nl]%*%solve(V[(i+1):nl,(i+1):nl])
		V_trans[i] <- V_trans[i]-tmp%*%V[(i+1):nl,i]
		theta_trans[i,] <- theta_hat[i,]+tmp%*%(theta[(i+1):nl,]-theta_hat[(i+1):nl,])
		S <- S+(theta[i,]-theta_trans[i,])%*%t(theta[i,]-theta_trans[i,])/V_trans[i]
	}
	}
	S <- S+(theta[nl,]-theta_trans[nl,])%*%t(theta[nl,]-theta_trans[nl,])/V_trans[nl]
	S/K
}

prcal <- function (theta_hat0,theta_hat,theta,omega) {
	nl <- length(theta_hat0)
	K <- dim(theta)[2]
	V <- matrix(NA,nl,nl)
	diag(V) <- theta_hat0*(1-theta_hat0)
	for (i in 1:(nl-1)) {
		for (j in (i+1):nl) {
			V[i,j] <- V[j,i] <- -theta_hat0[i]*theta_hat0[j]
		}
	}
	theta_trans <- theta_hat
	V_trans <- diag(V)
	if (nl>2) {
	for (i in 2:(nl-1)) {
		tmp <- V[i,(i+1):nl]%*%solve(V[(i+1):nl,(i+1):nl])
		V_trans[i] <- V_trans[i]-tmp%*%V[(i+1):nl,i]
		theta_trans[i,] <- theta_hat[i,]+tmp%*%(theta[(i+1):nl,]-theta_hat[(i+1):nl,])
	}
	}
	lnPr <- 0
	for (i in 2:nl) {
		lnPr <- lnPr + dmvnorm(x=theta[i,],mean=theta_trans[i,],sigma=V_trans[i]*omega,log=T)
	}
	lnPr
}

likcal <- function (Y,theta) {
    lnL <- 0
    if (max(Y[!is.na(Y)])>1) {
    for (i in 1:dim(Y)[1]) {
        if (!is.na(sum(Y[i,]))) {
            tmp <- dmultinom(Y[i,],size=sum(Y[i,]),prob=theta[,i],log=T)
            if (is.infinite(tmp) || is.na(tmp)) {
                tmp <- -1000
            }
            lnL <- lnL+tmp
        }
    }
    } else {
        tmp <- t(theta)*(Y==0)
        lnL <- sum(log(t(theta)[which(Y==1,arr.ind=T)]))+sum(rowSums(tmp==0)*log(1-rowSums(tmp)))
		if (is.infinite(lnL) || is.na(lnL))
		lnL <- -1000
    }
    lnL
}
