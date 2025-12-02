#update the WF solution to PIM process 09/09/2024

##Parameters to estimate: 
#m = number of memory updates. Its prior is a gamma distribution with mean=mu_m, so alpha=beta*mu_m, beta=rate_m.
#r = a vector with each element models the bias to use a variant type. If estimate.r=F, bias is not modelled. Prior for each r is an exponential distribution with rate_r.
#theta = a list, with each cell corresponding to a variable. Each cell contains a matrix of usage frequency of each variant by each speaker. Rows are variants and columns are speakers.
#omega = covariance matrix among speakers.
#beta = a vector of coefficients of each linguistic factor and each social factor, descrbing the amount of increase in the relative frequency of using the variant to a reference variant. The prior for each beta is normal with mean=0, sd=sd_beta. 

##Required Inputs:
#data
#data$Xind a matrix records social factors (row) of each speaker (col)
#data$Xinter a list, with each cell corresponding to to a variable. Each cell contains a matrix recording linguistic factors (row) of each speaker (col)
#data$lingua is a list, with each cell corresponding to a variable. Each cell contains a matrix, in which rows are speakers and columns are variants. Binary data is recorded as 0 and 1 for not use and use. Frequency data is recorded as integers, the number of times a variant is used.
#data$type is a list, with each cell corresponding to a variable. Each cell contains a vector, in which each element is the type of each variant in the same order as the variants in data$lingua. Heritage type = 1; Unrecorded type = 4. Unrecorded type is added for each variable with only one variant, such as comprehension variable.
#data$theta0 is a list, with each cell corresponding to a variable. Each cell contains a vector, in whcih each element is the initial usage frequency of each variant in the same order as the variants in data$lingua.
#inter = number of MCMC interations
#interval = number of intervals to write MCMCs to a file
#pathname = the path of the folder to save results
#mu_m = mean of gamma prior for m, also set as the starting value of m
#estimate.r = if bias is modelled
#eta_r = #tuning parameter for r
#eta_m = #tuning parameter for m
#eta_beta = #tuning parameter for each beta
#rate_r = #beta of gamma prior for each r, this give a very large variance to the prior
#rate_m = #beta of gamma prior for m, this give a very large variance to the prior
#sd_beta = #s.d. of normal prior for each beta, this give a very large variance to the prior

##Outputs are in saved files:
#post.csv = each row is a sampled MCMC after each interval, recording interation #,log-likelihood,log-prior,beta,r,m
#thetai.csv = the MCMC sample of theta for the ith language variable. Each sample is a matrix, with rows are variants and columns are speakers.
#omega.csv = the MCMC sample of the covariance matrix among speakers. Each sample is a KbyK matrix, where K is the number of speakers.

BayesVarbrul <- function (data,inter,interval,pathname,mu_m=100,estimate.r=FALSE,eta_r=0.001,eta_m=0.001,eta_beta=0.01,rate_r=1,rate_m=0.01,sd_beta=0.5) {
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
	K <- dim(data$Xind)[2] #K number of speakers
	Pind <- dim(data$Xind)[1] #Pind number of social factors for speakers
	#Pvar <- dim(data$Xindvar)[1] #Pvar number of factor for variants
	Pinter <- length(data$Xinter) #Pinter number of interactions
	n <- max(unlist(data$type)) #n number of variant types
	nl <- sapply(data$type,function (i) length(i)) #nl number of variants for each VARIABLE
	rho <- K+sum(nl)
	N <- 100 #memory size, default as a fixed value for model identifiability
	r <- rep(1/N,n)
	m <- mu_m
	R <- K*diag(1,K)
	lnL_list <- numeric(L) #list for log-likelihood, with each element corresponding to each language variable
	lnPr_list <- numeric(L) #llist for log-prior, with each element corresponding to each language variable
	theta_hat0 <- vector("list",L) #usage frequency expected from null model
	theta_hat <- vector("list",L) #usage frequency expected from null model and effects of social factors
	theta <- vector("list",L) 
	V_hat <- vector("list",L) 
	omega <- riwish(rho,R)
	S <- matrix(0,K,K)
	P <- Pind+Pinter
	beta <- numeric(P)
	for (i in 1:L) {
		r_tmp <- r[data$type[[i]]]
		for (z in 1:n) {
			idx <- which(data$type[[i]]==z)
			if (length(idx)>1) {
				r_tmp[idx] <- r_tmp[idx]/length(idx)
			}
		} 
		r0 <- sum(r_tmp)
		# introduce new mutation in the starting value
		theta_hat0[[i]] <- data$theta0[[i]]
		theta_hat0[[i]][theta_hat0[[i]]==1] <- 1-sum(theta_hat0[[i]]==0)/N
		theta_hat0[[i]][theta_hat0[[i]]==0] <- 1/N
		# use PIM model
		V_hat[[i]] <- 1/(1+2*r0*N)*(1-exp(-(1+2*r0*N)*m/N))*(diag(r_tmp/r0)-(r_tmp/r0)%*%t(r_tmp/r0))-(1-exp(-m/N))*exp(-(1+2*r0*N)*m/N)*(theta_hat0[[i]]-r_tmp/r0)%*%t(theta_hat0[[i]]-r_tmp/r0)+exp(-r0*m)/(1+r0*N)*(1+exp(-(1+r0*N)*m/N))*(diag(theta_hat0[[i]]-r_tmp/r0)-(theta_hat0[[i]]-r_tmp/r0)%*%t(r_tmp/r0)-(r_tmp/r0)%*%t(theta_hat0[[i]]-r_tmp/r0))
		theta_hat0[[i]] <- theta_hat0[[i]] %*%(exp(-r0*m)*(diag(nl[i])-rep(1,nl[i])%*%t(r_tmp/r0))+rep(1,nl[i])%*%t(r_tmp/r0))
		theta_hat[[i]] <- matrix(rep(theta_hat0[[i]],K),ncol=K,byrow=F)
		anc <- which(data$type[[i]]==1) #anc tells which variants are heritage type.
		if (Pinter>0) {
			for (z in 1:Pinter) {
				theta_hat[[i]] <- exp(log(theta_hat[[i]])+beta[Pind+z]*data$Xinter[[z]][[i]])
			}
			theta_hat[[i]] <- theta_hat[[i]]/t(matrix(rep(colSums(theta_hat[[i]]),nl[i]),K,nl[i]))
		}
		if (Pind>0) {
			if (nl[i]>length(anc)) {
				theta_hat[[i]][-anc,] <- theta_hat[[i]][-anc,]/sum(theta_hat[[i]][-anc])*matrix(rep(exp(log(sum(theta_hat[[i]][-anc])/sum(theta_hat[[i]][anc]))+beta[1:Pind]%*%data$Xind),nl[i]-length(anc)),nrow=nl[i]-length(anc),byrow=T)*sum(theta_hat[[i]][anc])
				theta_hat[[i]] <- theta_hat[[i]]/t(matrix(rep(colSums(theta_hat[[i]]),nl[i]),K,nl[i]))
			}
		}
		theta[[i]] <- theta_hat[[i]]
		lnL_list[i] <- likcal(data$lingua[[i]],theta[[i]])
		lnPr_list[i] <- prcal(theta_hat[[i]],V_hat[[i]],theta[[i]],omega)
	}
	lnL <- sum(lnL_list) #starting value of log likelihood
	
	if (estimate.r) {
		lnPr <- sum(lnPr_list)+sum(dexp(r,rate=rate_r,log=T))+sum(dnorm(beta,0,sd_beta,log=T))+dgamma(m,shape=rate_m,rate=rate_m,log=T) #starting value of log prior
	} else {
		lnPr <- sum(lnPr_list)+sum(dexp(r[1],rate=rate_r,log=T))+sum(dnorm(beta,0,sd_beta,log=T))+dgamma(m,shape=rate_m,rate=rate_m,log=T)
	}
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
			lnPr_tmp <- prcal(theta_hat[[i]],V_hat[[i]],theta_new,omega)
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
				if (Pinter>0) {
					for (z in 1:Pinter) {
						theta_hat_new[[i]] <- exp(log(theta_hat_new[[i]])+beta_new[Pind+z]*data$Xinter[[z]][[i]])
					}
					theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
				}
				if (Pind>0) {
					if (nl[i]>length(anc)) {
						theta_hat_new[[i]][-anc,] <- theta_hat_new[[i]][-anc,]/sum(theta_hat_new[[i]][-anc])*matrix(rep(exp(log(sum(theta_hat_new[[i]][-anc])/sum(theta_hat_new[[i]][anc]))+beta_new[1:Pind]%*%data$Xind),nl[i]-length(anc)),nrow=nl[i]-length(anc),byrow=T)*sum(theta_hat_new[[i]][anc])
						theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
					}
				}
				lnPr_list_new[i] <- prcal(theta_hat_new[[i]],V_hat[[i]],theta[[i]],omega)
				S_new <- S_new+scal(theta_hat_new[[i]],V_hat[[i]],theta[[i]])
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
		if (estimate.r) {
			#modeling bias, then updating each r
			for (j in 1:n) {
			u <- runif(1)
			r_new <- r
			r_new[j] <- r[j]*exp(eta_r*(u-0.5))
			lnPr_list_new <- numeric(L)
			S_new <- matrix(0,K,K)
			theta_hat_new <- vector("list",L)
			V_hat_new <- vector("list",L)
			theta_hat0_new <- vector("list",L)
			for (i in 1:L) {
				r_tmp <- r_new[data$type[[i]]]
				for (z in 1:n) {
				idx <- which(data$type[[i]]==z)
				if (length(idx)>1) {
					r_tmp[idx] <- r_tmp[idx]/length(idx)
				}
				} 
				r0 <- sum(r_tmp)
				# introduce new mutation in the starting value
				theta_hat0_new[[i]] <- data$theta0[[i]]
				theta_hat0_new[[i]][theta_hat0_new[[i]]==1] <- 1-sum(theta_hat0_new[[i]]==0)/N
				theta_hat0_new[[i]][theta_hat0_new[[i]]==0] <- 1/N
				# use PIM model
				V_hat_new[[i]] <- 1/(1+2*r0*N)*(1-exp(-(1+2*r0*N)*m/N))*(diag(r_tmp/r0)-(r_tmp/r0)%*%t(r_tmp/r0))-(1-exp(-m/N))*exp(-(1+2*r0*N)*m/N)*(theta_hat0_new[[i]]-r_tmp/r0)%*%t(theta_hat0_new[[i]]-r_tmp/r0)+exp(-r0*m)/(1+r0*N)*(1+exp(-(1+r0*N)*m/N))*(diag(theta_hat0_new[[i]]-r_tmp/r0)-(theta_hat0_new[[i]]-r_tmp/r0)%*%t(r_tmp/r0)-(r_tmp/r0)%*%t(theta_hat0_new[[i]]-r_tmp/r0))
				theta_hat0_new[[i]] <- theta_hat0_new[[i]] %*%(exp(-r0*m)*(diag(nl[i])-rep(1,nl[i])%*%t(r_tmp/r0))+rep(1,nl[i])%*%t(r_tmp/r0))

				theta_hat_new[[i]] <- matrix(rep(theta_hat0_new[[i]],K),ncol=K,byrow=F)
				anc <- which(data$type[[i]]==1)
				if (Pinter>0) {
					for (z in 1:Pinter) {
						theta_hat_new[[i]] <- exp(log(theta_hat_new[[i]])+beta[Pind+z]*data$Xinter[[z]][[i]])
					}
					theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
				}
				if (Pind>0) {
					if (nl[i]>length(anc)) {
						theta_hat_new[[i]][-anc,] <- theta_hat_new[[i]][-anc,]/sum(theta_hat_new[[i]][-anc])*matrix(rep(exp(log(sum(theta_hat_new[[i]][-anc])/sum(theta_hat_new[[i]][anc]))+beta[1:Pind]%*%data$Xind),nl[i]-length(anc)),nrow=nl[i]-length(anc),byrow=T)*sum(theta_hat_new[[i]][anc])
						theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
					}
				}
				lnPr_list_new[i] <- prcal(theta_hat_new[[i]],V_hat_new[[i]],theta[[i]],omega)
				S_new <- S_new+scal(theta_hat_new[[i]],V_hat_new[[i]],theta[[i]])
			}
			lnPr_new <- lnPr-sum(lnPr_list)+sum(lnPr_list_new)-dexp(r[j],rate=rate_r,log=T)+dexp(r_new[j],rate=rate_r,log=T)
			accept <- min(1,exp(lnPr_new-lnPr) * exp(eta_r*(u-0.5)))
			if (runif(1)<accept) {
				lnPr <- lnPr_new
				lnPr_list <- lnPr_list_new
				theta_hat[[i]] <- theta_hat_new[[i]]
				theta_hat0[[i]] <- theta_hat0_new[[i]]
				V_hat[[i]] <- V_hat_new[[i]]
				r <- r_new
				S <- S_new
			}
		}
		} else {
		#no bias, then all rs are the same
		u <- runif(1)
		r_new <- r*exp(eta_r*(u-0.5))
		lnPr_list_new <- numeric(L)
		S_new <- matrix(0,K,K)
		theta_hat_new <- vector("list",L)
		V_hat_new <- vector("list",L)
		theta_hat0_new <- vector("list",L)
		for (i in 1:L) {
			r_tmp <- r_new[data$type[[i]]]
			for (z in 1:n) {
			idx <- which(data$type[[i]]==z)
			if (length(idx)>1) {
				r_tmp[idx] <- r_tmp[idx]/length(idx)
			}
			} 
			
			r0 <- sum(r_tmp)
			# introduce new mutation in the starting value
			theta_hat0_new[[i]] <- data$theta0[[i]]
			theta_hat0_new[[i]][theta_hat0_new[[i]]==1] <- 1-sum(theta_hat0_new[[i]]==0)/N
			theta_hat0_new[[i]][theta_hat0_new[[i]]==0] <- 1/N
			# use PIM model
			V_hat_new[[i]] <- 1/(1+2*r0*N)*(1-exp(-(1+2*r0*N)*m/N))*(diag(r_tmp/r0)-(r_tmp/r0)%*%t(r_tmp/r0))-(1-exp(-m/N))*exp(-(1+2*r0*N)*m/N)*(theta_hat0_new[[i]]-r_tmp/r0)%*%t(theta_hat0_new[[i]]-r_tmp/r0)+exp(-r0*m)/(1+r0*N)*(1+exp(-(1+r0*N)*m/N))*(diag(theta_hat0_new[[i]]-r_tmp/r0)-(theta_hat0_new[[i]]-r_tmp/r0)%*%t(r_tmp/r0)-(r_tmp/r0)%*%t(theta_hat0_new[[i]]-r_tmp/r0))
			theta_hat0_new[[i]] <- theta_hat0_new[[i]] %*%(exp(-r0*m)*(diag(nl[i])-rep(1,nl[i])%*%t(r_tmp/r0))+rep(1,nl[i])%*%t(r_tmp/r0))

			theta_hat_new[[i]] <- matrix(rep(theta_hat0_new[[i]],K),ncol=K,byrow=F)
			anc <- which(data$type[[i]]==1)
			if (Pinter>0) {
				for (z in 1:Pinter) {
					theta_hat_new[[i]] <- exp(log(theta_hat_new[[i]])+beta[Pind+z]*data$Xinter[[z]][[i]])
				}
				theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
			}
			if (Pind>0) {
				if (nl[i]>length(anc)) {
					theta_hat_new[[i]][-anc,] <- theta_hat_new[[i]][-anc,]/sum(theta_hat_new[[i]][-anc])*matrix(rep(exp(log(sum(theta_hat_new[[i]][-anc])/sum(theta_hat_new[[i]][anc]))+beta[1:Pind]%*%data$Xind),nl[i]-length(anc)),nrow=nl[i]-length(anc),byrow=T)*sum(theta_hat_new[[i]][anc])
					theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
				}
			}
			lnPr_list_new[i] <- prcal(theta_hat_new[[i]],V_hat_new[[i]],theta[[i]],omega)
			S_new <- S_new+scal(theta_hat_new[[i]],V_hat_new[[i]],theta[[i]])
		}
		lnPr_new <- lnPr-sum(lnPr_list)+sum(lnPr_list_new)-dexp(r[1],rate=rate_r,log=T)+dexp(r_new[1],rate=rate_r,log=T)
		accept <- min(1,exp(lnPr_new-lnPr)*exp(eta_r*(u-0.5)))
		if (runif(1)<accept) {
			lnPr <- lnPr_new
			lnPr_list <- lnPr_list_new
			theta_hat[[i]] <- theta_hat_new[[i]]
			theta_hat0[[i]] <- theta_hat0_new[[i]]
			V_hat[[i]] <- V_hat_new[[i]]
			r <- r_new
			S <- S_new
		}
		}
		#updating m
		u <- runif(1)
		m_new <- m*exp(eta_m*(u-0.5))
		lnPr_list_new <- numeric(L)
		S_new <- matrix(0,K,K)
		theta_hat_new <- vector("list",L)
		V_hat_new <- vector("list",L)
		theta_hat0_new <- vector("list",L)
		for (i in 1:L) {
			r_tmp <- r[data$type[[i]]]
			for (z in 1:n) {
				idx <- which(data$type[[i]]==z)
				if (length(idx)>1) {
					r_tmp[idx] <- r_tmp[idx]/length(idx)
				}
			}
			r0 <- sum(r_tmp)
			# introduce new mutation in the starting value
			theta_hat0_new[[i]] <- data$theta0[[i]]
			theta_hat0_new[[i]][theta_hat0_new[[i]]==1] <- 1-sum(theta_hat0_new[[i]]==0)/N
			theta_hat0_new[[i]][theta_hat0_new[[i]]==0] <- 1/N
			# use PIM model
			V_hat_new[[i]] <- 1/(1+2*r0*N)*(1-exp(-(1+2*r0*N)*m/N))*(diag(r_tmp/r0)-(r_tmp/r0)%*%t(r_tmp/r0))-(1-exp(-m/N))*exp(-(1+2*r0*N)*m/N)*(theta_hat0_new[[i]]-r_tmp/r0)%*%t(theta_hat0_new[[i]]-r_tmp/r0)+exp(-r0*m)/(1+r0*N)*(1+exp(-(1+r0*N)*m/N))*(diag(theta_hat0_new[[i]]-r_tmp/r0)-(theta_hat0_new[[i]]-r_tmp/r0)%*%t(r_tmp/r0)-(r_tmp/r0)%*%t(theta_hat0_new[[i]]-r_tmp/r0))
			theta_hat0_new[[i]] <- theta_hat0_new[[i]] %*%(exp(-r0*m)*(diag(nl[i])-rep(1,nl[i])%*%t(r_tmp/r0))+rep(1,nl[i])%*%t(r_tmp/r0))

			theta_hat_new[[i]] <- matrix(rep(theta_hat0_new[[i]],K),ncol=K,byrow=F)
			anc <- which(data$type[[i]]==1)
			if (Pinter>0) {
				for (z in 1:Pinter) {
					theta_hat_new[[i]] <- exp(log(theta_hat_new[[i]])+beta[Pind+z]*data$Xinter[[z]][[i]])
				}
				theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
			}
			if (Pind>0) {
				if (nl[i]>length(anc)) {
					theta_hat_new[[i]][-anc,] <- theta_hat_new[[i]][-anc,]/sum(theta_hat_new[[i]][-anc])*matrix(rep(exp(log(sum(theta_hat_new[[i]][-anc])/sum(theta_hat_new[[i]][anc]))+beta[1:Pind]%*%data$Xind),nl[i]-length(anc)),nrow=nl[i]-length(anc),byrow=T)*sum(theta_hat_new[[i]][anc])
					theta_hat_new[[i]] <- theta_hat_new[[i]]/t(matrix(rep(colSums(theta_hat_new[[i]]),nl[i]),K,nl[i]))
				}
			}
			lnPr_list_new[i] <- prcal(theta_hat_new[[i]],V_hat_new[[i]],theta[[i]],omega)
			S_new <- S_new+scal(theta_hat_new[[i]],V_hat_new[[i]],theta[[i]])
		}
		lnPr_new <- lnPr-sum(lnPr_list)+sum(lnPr_list_new)-dgamma(m,shape=rate_m,rate=rate_m,log=T)+dgamma(m_new,shape=rate_m,rate=rate_m,log=T)
		accept <- min(1,exp(lnPr_new-lnPr)*exp(eta_m*(u-0.5)))
		if (runif(1)<accept) {
			lnPr <- lnPr_new
			lnPr_list <- lnPr_list_new
			theta_hat[[i]] <- theta_hat_new[[i]]
			theta_hat0[[i]] <- theta_hat0_new[[i]]
			V_hat[[i]] <- V_hat_new[[i]]
			m <- m_new
			S <- S_new
		}
		#updating omega
		omega <- riwish(rho,R+S)
		for (i in 1:L) {
			lnPr_list[i] <- prcal(theta_hat[[i]],V_hat[[i]],theta[[i]],omega)
		}
		if (estimate.r) {
		lnPr <- sum(lnPr_list)+sum(dexp(r,rate=rate_r,log=T))+sum(dnorm(beta,0,sd_beta,log=T))+dgamma(m,shape=rate_m,rate=rate_m,log=T) #starting value of log prior
	} else {
		lnPr <- sum(lnPr_list)+sum(dexp(r[1],rate=rate_r,log=T))+sum(dnorm(beta,0,sd_beta,log=T))+dgamma(m,shape=rate_m,rate=rate_m,log=T)
	}
			
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
scal <- function (theta_hat,V_hat,theta) {
	nl <- dim(theta)[1]
	K <- dim(theta)[2]
	V <- V_hat
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

prcal <- function (theta_hat,V_hat,theta,omega) {
	nl <- dim(theta)[1]
	K <- dim(theta)[2]
	V <- V_hat
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
