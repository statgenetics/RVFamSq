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
    
    
    if(ncol(simdata$covar)>1) {
      EYi=miu+rowSums(t(apply(matrix(simdata$covar[sid,],length(sid),length(betax)),1,function(x) x*betax)))
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
Tscore<-function(ped,maf,maf_cutoff,paras,covar_update_col,trait_update_col,pop_update_col,geno_update_col) {
  miu<-paras$miu
  sitag<-paras$sitag
  sitae<-paras$sitae
  betax<-unlist(paras[which(names(paras)!="miu"&names(paras)!="sitag"&names(paras)!="sitae")])
  
  if (is.null(nrow(ped))) 
    ped<-matrix(ped,1,length(ped))
  
  covs<-as.numeric(ped[,covar_update_col])
  covs<-matrix(covs,nrow(ped),length(covar_update_col))
  
  if(!is.null(ncol(ped[,covar_update_col]))) {
    EYi=miu+rowSums(t(apply(covs,1,function(x) x*betax)))
  } else {
    EYi=miu+betax*as.numeric(ped[,covar_update_col])
  }
  
  maf_pop<-as.matrix(t(cbind(maf[,as.numeric(ped[,pop_update_col])])))
  RVs_cutoff<-which(maf_pop>maf_cutoff,arr.ind = TRUE)
  
  
  yi<-as.numeric(ped[,trait_update_col])
  kindex<- match(ped[,2],rownames(kin))
  fanij<-as.matrix(2*kin[kindex,kindex]*sitag^2)
  
  diag(fanij)<-sitag^2 + sitae^2
  fanij_inv<-solve(fanij)
  
  genos<-matrix(as.numeric(ped[,geno_update_col]),nrow=nrow(ped),ncol = length(geno_update_col))
  genos[RVs_cutoff]<-0
  
  g<-rowSums(genos) - rowSums(2*maf_pop)
  
  upperdot = t(g) %*% fanij_inv %*% (yi - EYi)
  lowerdot = t(g) %*% fanij_inv %*% g
  return(c(upperdot,lowerdot))
}
