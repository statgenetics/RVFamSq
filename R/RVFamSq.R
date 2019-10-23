#' Rare Variant-Family-Based Score Test for Quantitative Traits
#' 
#' A regional association analysis of rare variants and quantitative traits in
#' family data
#'RV_FamSq <- function(ped_pheno, ped_geno, maf_data, maf_cutoff,parafile,covar_col, trait_col, pop_col=c(), out,kin, estimateAF,start_par)
#' @param ped_pheno a dataframe of phenotype. The first four columns should be family IDs, 
#' individual IDs, father IDs and mother IDs of each sample.
#' The ped_pheno can also include columns of covariates and populations IDs. 
#' The populations IDs indicate which corresponding MAF will be used in the model. 
#' @param ped_geno a dataframe of genotype. The first two columns are family IDs and individual IDs. 
#' The rest columns of the ped_geno are genotype score at each sites. The value of genotype score is 0, 1 or 2,
#' corresponding to number of rare variants observed at each sites. Number of sites included in ped_geno should be equal to the number of rows 
#' of \code{maf_data}. 
#' @param maf_data a dataframe describing the information of rare variants included in ped_geno.
#' The first three columns correspond to the name of the gene, chromosome number, and position of variants. The rest columns of maf_data
#' are MAFs of variants in each specific population group obtained from database, e.g., Genome Aggregation Database (gnomAD).
#' @param maf_cutoff the cutoff of MAF for rare variants 
#' @param parafile the path to the file saving parameters estimated under the null model.
#' If the file is not existed, the package will estimate the parameters
#' using the available data and create a parafile in the path.
#' @param covar_col a integer or a vector of integers indicating which columns of \code{ped_pheno} will be used as covariants of the analysis. Note that 
#' multiple covariants are allowed in the analysis. 
#' @param trait_col a integer indicating which columns of \code{ped_pheno} will be used as traits.
#' @param pop_col a integer indicating which columns of \code{ped_pheno} describe the population group of the sample.
#' @param out the directory that save the results. RVFamSq outputs a dataframe results including
#' the name of the gene, score, p-value of the gene, number of sample size
#' and number of families
#' @param kin a kinship matrix calculated either based on family strtuctures or genomic informations. 
#' @param estimateAF a dataframe of estimated MAFs that do not oberserved in public datasets. The headers decribing the name 
#' of the sub-population should be consistent with the headers in \code{maf_data}. If this argument is not defined, RVfamSq will
#' estimate the missing MAF based on the dosage of founder. 
#' @param start_par an optional argument that defines the starting values of
#' the parameters fitting. User can define specific values to start the
#' parameters fitting, otherwise, the starting values will generate randomly.
#' @return \item{results}{a data frame results containing the name of the gene,
#' the value of the score and p-value write into the file named by the name of
#' the gene under the directory defined by the argument \code{out}.}
#' \item{paras}{a list of parameters estimated by maximizing the likelihood.
#' The values of the paras save to the file defined by the argument
#' \code{parafile}. If the \code{parafile} is existed, the values of the paras
#' load from the file \code{parafile}.}
#' @author Zhihui Zhang and Suzanne Leal
#' @details RVFamSq performed association analysis by estimating the parameters under the null model and calculating the statistical score using these parameters.
#' The parameters are estimated by maximizing the multivariate normal likelihood: 
#' \deqn{L=\prod_{i}(2\pi)^{-n_i / 2}\left |\Omega _i\right |^{-1 /2}e^{\left [ y_i - E(y_i)\right ]'\Omega _i^{-1}\left [ y_i - E(y_i)\right ]}} 
#' where \eqn{n_i} is the number of individuals in family i and \eqn{\left |\Omega _i\right |} is the determinant of matrix \eqn{\Omega _i}. \eqn{E(y_i)} in the above
#' function is a vector of expected phenotype of all individuals within each family and is calculated by: 
#' \deqn{E(y_i)=\mu + \beta _x x_i}
#' where \eqn{x_i} is a \eqn{n_i \times q} matrix of q covariates included in the association analysis. The parameters \eqn{\mu} and \eqn{\beta _x} are population mean and a vector of covariate effects. 
#' The two parameters are estimated by maximizing the above likelihood function. The matrix \eqn{\Omega_i} in the likelihood function is calculated within each family i and defined as:
#' \deqn{\Omega  _{ijk}\left\{\begin{matrix}\sigma _g^2 + \sigma_e^2\ \ \ \ \ \ \ \ \ \ \ \ \ if\  j=k\\
#' 2\varphi _{ijk}\sigma _g^2\ \ \ \ \ \ \ \ \ \ \ \ if\ i\neq k 
#' \end{matrix}\right.}
#' Here, the parameter \eqn{\sigma_g^2} and\eqn{\sigma_e^2} are variance components that are account for background polygenic effects
#' and environmental effects, respectively. Additionally, \eqn{\Omega _{ijk}} denote the kinship coefficient between individuals j and k. 
#' The values of parameters \eqn{\sigma_g^2}, \eqn{\sigma_e^2}, \eqn{\beta_x}, and \eqn{\mu} are estimated by maximizing the above likelihood.
#' With the values of these parameters, we can obtain the variance-covariance matrix  \eqn{\Omega _i^{base}} and vector \eqn{E(y_i)^{base}} based on the above equations for each family. 
#' Using these two quantities, the statistic score is defined by:
#' \deqn{T^{SCORE}=\frac{\left \{ \sum _i \left [G_i -E(G_i)\right ]'[\Omega _i^{(base)}]^{-1}[y_i-E(y_i)^{(base)}] \right \}^2}{\sum _i[G_i-E(G_i)]'[\Omega_i^{(base)}]^{-1}[G_i-E(G_i)^{(base)}] }}
#' where \eqn{E(G_i)} is defined as:
#' \deqn{E(G_i)=\sum_{m=1}^{M}2p_m}
#' where \eqn{p_m} is the minor allele frequency (MAF) which is obtain from the dataframe of \code{maf_data}. 
#' @section Data format:
#' Example data is generated for 1,000 extended families with 10 members in each family. The information of RVs within functional region of the gene AAAS with MAF<2% are obtained from gnomeAD. 
#' Quantitative traits are generated randomly based on the disribution of N(0,1) and genotypes of samples are generated based on MAF of the gene AAAS.
#' \itemize{
#' \item{\strong{ped_pheno.txt}} {A phenotype file of 1,000 extended families (10,000 individuals in total). Each column of the file represents family IDs, individual IDs, father IDs, mother IDs, sex, traits, and population IDs, respectively. 
#' RVFamSq is capable to analyze multiple covariates and extra covariants should be also included in this file after the fourth column.}
#' \item{\strong{ped_geno.txt}} {A genotype file of samples included in the ped_pheno.txt file. The first two columns of the file should be family IDs and individual IDs that have the same format as in ped_pheno.txt. 
#' The rest columns of the file are genotypes of each individual and codes as 0, 1, and 2.}
#' \item{\strong{AAAS.sfs}} {A file with descriptive information of RVs within the functional region of the gene AAAS. The file contains four columns representing the name of the interested gene, chromosome, position, and MAF of the variant, respectively.
#' If more than one sub-populations are included in the samples, multiple MAFs correspond to each specific population should be included in this file and denoted by headers of population specific symbols}
#' }
#' @references Chen, W.-M., and Abecasis, G.R. (2007). Family-Based Association
#' Tests for Genomewide Association Scans. Am. J. Hum. Genet. 81, 913–926.\cr
#' @examples
#' 
#' library (kinship2)
#' library (RVFamSq)
#' 
#' ##  load example data
#' data_dir<-system.file("data", package = "RVFamSq")
#' ped_pheno<-read.table(phenofile<-paste0(data_dir,"/ped_pheno.txt"))
#' ped_geno<-read.table(paste0(data_dir, "/ped_geno.txt"))
#' maf_data<-read.table(paste0(data_dir,"/AAAS.sfs"), header = TRUE)
#' 
#' ## Define the files that save the parameters estimated under the null model.
#' ## If the file is not existed, the package will estimate the parameters based on the available data.
#' parafile<-paste0(data_dir,"/results/paras.rds")
#' 
#' ##  Define output directory that save the results
#' out<-paste0(data_dir,"/results")
#' 
#' ## Load phenotype data and calculate the kinship matrix of the pedigree. In this example, we utilize the R package to estimate the kinship of samples based on family structure, but the kinship can also be estimated by genetic variants.
#' kin_pre<-data.frame(id=ped_pheno[,2], mom=ped_pheno[,4], dad=ped_pheno[,3], sex=ped_pheno[,5])
#' tped <- with(kin_pre, pedigree(id, dad, mom, sex, famid=ped_pheno[,1]))
#' kin <- kinship(tped)
#' 
#' ## Run RVFamsSq package and calculate the statistical score of the interested gene.
#' RV_FamSq(ped_pheno=ped_pheno, ped_geno=ped_geno, maf_data=maf_data, maf_cutoff =0.02,parafile=parafile,covar_col=c(5), trait_col=c(6), pop_col=c(7),out=out, kin=kin)
#'
#' @import bbmle 
#' @import mvtnorm 
#' @import rlist
#' @export
RV_FamSq <- function(ped_pheno, ped_geno, maf_data, maf_cutoff,parafile,covar_col, trait_col, pop_col=c(), out,kin, estimateAF,start_par) {
  library(bbmle)
  library(mvtnorm)
  library(rlist)
  
  #change dataframe to matrix by accelerating speed
  ped_pheno<-as.matrix(ped_pheno)
  ped_geno<-as.matrix(ped_geno)
  maf_data<-as.matrix(maf_data)
  
  if ((ncol(ped_geno)-2)!=nrow(maf_data)) stop("Dimension of genotype and number of varints in MAF file does not match")
  
  #choose RVs based on cutoff of MAF
  maf_pre<-matrix(as.numeric(maf_data[,4:ncol(maf_data)]),nrow=nrow(maf_data),ncol=ncol(maf_data)-3)
  rvs_index<-unique(which(maf_pre<maf_cutoff | is.na(maf_pre),arr.ind = TRUE)[,1])
  maf<-maf_pre[rvs_index,]
  colnames(maf)<-colnames(maf_data)[4:ncol(maf_data)]
  ped_geno<-ped_geno[,c(1,2,2+rvs_index)]
  
  #Delete samples without phenotypes or covariants
  miss_index<-which(is.na(ped_pheno[,covar_col]) | is.na(ped_pheno[,trait_col]))
  if(length(miss_index)>0) {
    print("Delete samples without the data of phenotype and covariant")
    ped_pheno<-ped_pheno[-miss_index,]
  }
  ped_pheno<-ped_pheno[,c(1:4,covar_col,trait_col,pop_col)]
  
  #Estimate mising MAF using external source or by dosage of founders
  if (length(which(is.na(maf))>0)) {
    na_index<-which(is.na(maf),arr.ind = TRUE)
    if (!missing(estimateAF)) {
      print("Estimate missing MAF using external 'estimateAF'")
      af_index<-sapply(colnames(maf), function(x) grep(x,colnames(estimateAF)))
      maf[na_index]<-as.numeric(estimateAF[af_index[na_index[,2]]])
    } else {
      print("Estimate missing MAF using founders' dosage.")
      founders<-ped_pheno[which((is.na(ped_pheno[,3]) & is.na(ped_pheno[,4])) | (ped_pheno[,3]==0 & ped_pheno[,4]==0)),2]
      founders_index<-match(founders, ped_geno[,2])
      na_ind<-which(!is.na(ped_geno[founders_index,3:ncol(ped_geno)]),arr.ind = TRUE)

      if (nrow(na_ind)>0) {
          print(paste(c("Estimate MAF based on dosage from", nrow(na_ind),"founders."),collapse = " "))
          maf[na_index]<-colSums(matrix(as.numeric(ped_geno[cbind(founders_index[na_ind[,1]],2+na_ind[,2])]),nrow=nrow(na_ind)))/(2*nrow(na_ind))
        } else {
          stop("No founders to estimate missing MAF, please refer optional methods to estimate missing MAF and input it as 'estimateAF'")
      }
    }
  }
  
  #Merge phenotype and genotype of the data
  ped_data<-as.matrix(merge(ped_pheno,ped_geno,by.x=c(1,2),by.y=c(1,2)))
    
  #get updated column number for the new merging ped_data
  covar_update_col<-c(5:(4+length(covar_col)))
  trait_update_col<-c(5+length(covar_col))
  pop_update_col<-c((6+length(covar_col)))
  geno_update_col<-c((7+length(covar_col)):ncol(ped_data))
  
  #assign geno score to missing genotypes based on binormial distribution using MAF
  if (length(which(is.na(ped_data[,geno_update_col]))>0)) {
    na_geno_index<-which(is.na(ped_data[,geno_update_col]),arr.ind = TRUE)
    maf_list<-maf[cbind(na_geno_index[,2],as.numeric(ped_data[na_geno_index[,1],pop_update_col]))]
    ped_data[cbind(na_geno_index[,1],geno_update_col[na_geno_index[,2]])]<-sapply(maf_list, function(x) sum(rbinom(2,1,x)))
  }
  
  genos<-matrix(as.numeric(ped_data[,geno_update_col]),ncol=nrow(maf),nrow = nrow(ped_data))
  
  if(missing(kin)) {
    stop("Kinship matrix is missing.")
  } else {
    if (rownames(kin) %contain% ped_data[,2]) {
      assign("kin", kin, envir=.GlobalEnv)
    } else {
      stop("Some samples in pedfile are not included in Kinship matrix")
    }
  }
  
  
  if(!file.exists(parafile)) {
    par_num<-3+length(covar_col)
    if(missing(start_par)) {
      print("Starting values of parameters are not defined, generating random initial values now...")
      start_par<-runif(par_num,-1,1)
    } else {
      if (length(start_par)<par_num) {
        print("Defined less parameters than required by the model, generating random initial values now...")
        start_par<-runif(par_num,-1,1)
      }
    }
    
    if (length(covar_col)>1) {
      assign("betax_names",paste0("betax",1:length(covar_col)), envir=.GlobalEnv)
    } else {
      assign("betax_names","betax", envir=.GlobalEnv)
    }
    
    paras <- start_par
    names(paras)<-c("miu",betax_names,"sitag","sitae")
    
    covs<-as.numeric(ped_data[,covar_update_col])
    covs<-matrix(covs,nrow(ped_data),length(covar_update_col))
    
    assign("simdata", list(id=ped_data[,2],pheno=as.numeric(ped_data[,trait_update_col]),covar=covs,famids=ped_data[,1],genos=as.numeric(genos), maf=as.numeric(maf)), envir=.GlobalEnv)
    print(paste("Estimating the values of parameters",paste(names(paras),collapse = ","),sep = " "))
    nullmod<-nlf_fit(paras)
    paras<-as.list(coef(nullmod))
    list.save(paras, parafile)
  } else {
    print("Reading parameters from parafile... ")
    paras<-list.load(parafile)
    print(unlist(paras))
    
    peds_grp_by_fam<-split(seq_len(nrow(ped_data)),ped_data[,1])
    scores<-unlist(lapply(peds_grp_by_fam, function(x) Tscore(ped_data[unlist(x),],maf,maf_cutoff,paras,covar_update_col,trait_update_col,pop_update_col,geno_update_col)))
    upper<-sum(scores[seq(1,length(scores),2)])
    lower<-sum(scores[seq(2,length(scores),2)])
    score<-upper^2/lower
    p=pchisq(score,df=1)
    
    gene_name=as.character(maf_data[1,1])
    out<-paste(c(out,"/", gene_name,".out"),collapse ="")
    pout<-data.frame(gene=gene_name, score=score, p=1-p, sample_size=nrow(ped_data), family_size=length(unique(ped_data[,1])))
    write.table(pout,out, row.names = FALSE, col.names = TRUE, quote=FALSE,sep="\t")
    print(pout)
  }
  return()
}
