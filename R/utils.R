"%contain%" <- function(values,x) {
  tx <- table(x)
  tv <- table(values)
  z <- tv[names(tx)] - tx
  all(z >= 0 & !is.na(z))
}

#define negative likelihood function
nlf <- function(paras) {
  miu<-paras["miu"]
  sitag<-paras["sitag"]
  sitae<-paras["sitae"]

  betax<-paras[which(names(paras)!="miu"&names(paras)!="sitag"&names(paras)!="sitae")]

  sum_likelihood=0
  for (n in unique(simdata$famid)) {
    sid=which(simdata$famids==n)
    sid_count=length(sid)

    if(!is.null(ncol(simdata$covar))) {
      EYi=miu+rowSums(t(apply(simdata$covar[sid,],1,function(x) x*betax)))
    } else {
      EYi=miu+betax*simdata$covar[sid]
    }
    yi=simdata$pheno[sid]

    fid<-simdata$id[sid]
    kindex<-match(fid,rownames(kin))

    fani<-as.matrix(2*kin[kindex,kindex]*sitag^2)
    diag(fani)<-sitag^2 + sitae^2

    #negative log likelihood
    if (is.finite(log(det(fani)))) {
      sum_likelihood=sum_likelihood+dmvnorm(yi, EYi,fani,log=TRUE)
    }
  }
  return(-2*sum_likelihood)
}

#maximum likelihood function to get parameters
nlf_fit<-function(start) {
  parnames(nlf)<-c("miu",betax_names,"sitag","sitae")
  m1 <- mle2(nlf,start=start, data=simdata, method = "Nelder-Mead", skip.hessian = FALSE, control=list(trace=FALSE, maxit=1000))
  m2 <- mle2(nlf,start=coef(m1), control=list(trace=FALSE, maxit=1000, parscale=coef(m1)), data=simdata)
  return(m2)
}

#main function of the score model
Tscore<-function(ped,maf,paras,covar_update_col,trait_update_col,geno_update_col) {
  miu<-paras$miu
  sitag<-paras$sitag
  sitae<-paras$sitae
  betax<-unlist(paras[which(names(paras)!="miu"&names(paras)!="sitag"&names(paras)!="sitae")])

  if(!is.null(ncol(ped[,covar_update_col]))) {
    EYi=miu+rowSums(t(apply(ped[,covar_update_col],1,function(x) x*betax)))
  } else {
    EYi=miu+betax*ped[,covar_update_col]
  }

  yi<-ped[,trait_update_col]
  kindex<- match(ped[,2],rownames(kin))

  fanij<-as.matrix(2*kin[kindex,kindex]*sitag^2)
  diag(fanij)<-sitag^2 + sitae^2
  fanij_inv<-solve(fanij)

  g<-rowSums(ped[,geno_update_col]) - sum(2*maf)

  upperdot = t(g) %*% fanij_inv %*% (yi - EYi)
  lowerdot = t(g) %*% fanij_inv %*% g

  return(c(upperdot,lowerdot))
}
