sen <- function (y, z, mset, gamma = 1, inner = NULL, trim = NULL, lambda = 1/2,
              tau = 0, alternative="greater")
    {
        #Check input
        stopifnot((alternative=="greater")|(alternative=="less"))
        stopifnot(is.vector(gamma)&(length(gamma)==1))
        stopifnot(is.vector(lambda)&(length(lambda)==1))
        stopifnot(is.vector(tau)&(length(tau)==1))
        stopifnot(gamma>=1)
        stopifnot((lambda>0)&(lambda<1))
        stopifnot(is.vector(y)&is.vector(z)&is.vector(mset))
        stopifnot((length(z)==length(y)))
        stopifnot((length(z)==length(mset)))
        stopifnot(is.numeric(mset))
        stopifnot(all(!is.na(y)))
        stopifnot(all((z==0)|(z==1))) #z is 1 for treated, 0 for control
        tbcheck<-table(z,mset)
        ck<-all(tbcheck[2,]==1)&all(tbcheck[1,]>=1)
        if (!ck){
          warning("Every matched set must contain one treated subject and at least one control.")
          stopifnot(ck)
        }

        #Convert y to matrix
        mset<-as.integer(mset)
        o<-order(mset,1-z)
        y<-y[o]
        z<-z[o]
        mset<-mset[o]
        tb<-table(mset) #need to check
        nset<-length(tb)
        setsize<-max(tb)

        # Set default value of inner
        if (is.null(inner)){
          it<-c(.8,.8,.6,.4)
          cnts<-setsize-1
          if (cnts<5) inner<-it[1:cnts]
          else inner<-c(it,rep(0,cnts-4))
        }
        # Set default value of trim
        if (is.null(trim)) trim<-rep(3,setsize-1)

        # Vector expansion of scalar inner or trim
        if (is.vector(trim)&(length(trim)==1)) trim<-rep(trim,setsize-1)
        if (is.vector(inner)&(length(inner)==1)) inner<-rep(inner,setsize-1)

        # Check inner and trim
        stopifnot(is.vector(inner))
        stopifnot(is.vector(trim))
        stopifnot(length(trim)==length(inner))
        stopifnot(length(trim)==(setsize-1))
        stopifnot(all(inner>=0)&all(inner<=trim)&all(trim<Inf))

        makeymat<-function(yj){
          ymat<-matrix(NA,nset,setsize)
          m<-0
          for (i in 1:nset){
            ymat[i,1:tb[i]] <- yj[(m+1):(m+tb[i])]
            m<-m+tb[i]
          }
          ymat
        }

        ymat<-makeymat(y)

        if (alternative=="less"){
          ymat<-(-ymat)
          tau<-(-tau)
        }
        if (!(tau == 0)) ymat[, 1] <- ymat[, 1] - tau

        # sort ymat by numer of controls
        nct<-apply(!is.na(ymat),1,sum)-1 #number of controls
        o<-order(nct)
        ymat<-ymat[o,]
        nct<-nct[o]

        # Find scale factor
        n <- dim(ymat)[1]
        m <- dim(ymat)[2]
        out <- matrix(NA, n, m)
        one <- rep(1, m - 1)
        difs <- array(NA, c(n, m, m - 1))
        for (j in 1:m) {
          difs[, j, ] <- outer(as.vector(unlist(ymat[, j])), one,
                               "*") - ymat[, -j]
        }
        hqu <- as.numeric(stats::quantile(abs(as.vector(difs)), lambda
                                   , na.rm = TRUE))
        if (hqu<=0){
          warning("Error: Scale factor is zero.  Increase lambda.")
          stopifnot(hqu>0)
        }
        ymat<-ymat/hqu

        ms<-NULL
        u<-sort(unique(nct))
        for (i in 1:length(u)){
          j<-u[i]
          ms<-rbind(ms,mscoreInternal(ymat[nct==j,],inner[j],trim[j])/(j+1))
        }

        separable1v(ms, gamma = gamma)
    }
