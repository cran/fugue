senCI<-function(y, z, mset, gamma = 1, inner = NULL, trim = NULL, lambda = 1/2, 
    alpha = 0.05, alternative = "greater"){
	
	stopifnot((alternative=="greater")|(alternative=="less")|(alternative=="twosided"))
	stopifnot(is.vector(gamma) & (length(gamma) == 1))
	stopifnot(is.vector(lambda) & (length(lambda) == 1))
	stopifnot(is.vector(alpha) & (length(alpha) == 1))
	stopifnot(gamma >= 1)
	stopifnot((lambda > 0) & (lambda < 1))
	stopifnot((alpha>0) & (alpha<1))
	stopifnot(is.vector(y) & is.vector(z) & is.vector(mset))
	stopifnot((length(z) == length(y)))
	stopifnot((length(z) == length(mset)))
	stopifnot(is.numeric(mset))
	stopifnot(all(!is.na(y)))
	stopifnot(all((z == 0) | (z == 1)))
	tbcheck <- table(z, mset)
	ck <- all(tbcheck[2, ] == 1) & all(tbcheck[1, ] >= 1)
	if (!ck) {
	  warning("Every matched set must contain one treated subject and at least one control.")
	  stopifnot(ck)
    }
    if (alternative=="twosided") alpha<-alpha/2
	mx<-max(y)-min(y)
	int<-c(-mx,mx)
	crit<-qnorm(1-alpha)
	upper<-function(taus){
		nt<-length(taus)
		o<-rep(NA,nt)
		for (i in 1:nt){
			o[i]<-sen(y,z,mset,gamma=gamma,inner=inner,trim=trim,lambda=lambda,tau=taus[i],alternative="greater")$deviate-crit
			}
		o
	}
	lower<-function(taus){
		nt<-length(taus)
		o<-rep(NA,nt)
		for (i in 1:nt){
			o[i]<-sen(y,z,mset,gamma=gamma,inner=inner,trim=trim,lambda=lambda,tau=taus[i],alternative="less")$deviate-crit
			}
		o
	}
	uppere<-function(taus){
		nt<-length(taus)
		o<-rep(NA,nt)
		for (i in 1:nt){
			o[i]<-sen(y,z,mset,gamma=gamma,inner=inner,trim=trim,lambda=lambda,tau=taus[i],alternative="greater")$deviate
			}
		o
	}
	lowere<-function(taus){
		nt<-length(taus)
		o<-rep(NA,nt)
		for (i in 1:nt){
			o[i]<-sen(y,z,mset,gamma=gamma,inner=inner,trim=trim,lambda=lambda,tau=taus[i],alternative="less")$deviate
			}
		o
	}
	pe<-c(NA,NA)
	ci<-c(NA,NA)
	pe[1]<-stats::uniroot(uppere,int)$root
	pe[2]<-stats::uniroot(lowere,int)$root
	names(pe)<-c("lower","upper")
	names(ci)<-c("lower","upper")
	if (alternative=="greater"){
		ci[2]<-Inf
		ci[1]<-stats::uniroot(upper,int)$root
		}
	if (alternative=="less"){
		ci[1]<-(-Inf)
		ci[2]<-stats::uniroot(lower,int)$root
		}
	if (alternative=="twosided"){
		ci[1]<-stats::uniroot(upper,int)$root
		ci[2]<-stats::uniroot(lower,int)$root
		}
	list(point.estimates=pe,confidence.interval=ci)
}