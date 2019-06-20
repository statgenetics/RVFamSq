#' Rare Variant-Family-Based Score Test for Quantitative Traits
#' 
#' A regional association analysis of rare variants and quantitative traits in
#' family data
#'
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param genofile referring to the path of the genotype file. The columns in
#' the file represent family IDs, individual IDs and genotypes of each
#' individual. The file should not have header.
#' @param phenofile referring to the path of the phenotype file. The first five
#' columns in the file should be family IDs,individual IDs, father IDs, mother
#' IDs, and sex, respectively. The rest columns of the file are covariates and
#' trait. The covariates may be multi-columns and trait should be one column.
#' The file should not have header.
#' @param maffile the path to the .sfs file describing the information of
#' variants in the genotype file.
#' @param parafile the path to the file saving parameters estimated in null
#' model. If the file is not existed, the package will estimate the parameters
#' using the available data and save the values of the parameters in parafile.
#' @param covar_col a vector defining the index of columns of covariates in
#' \code{phenofile}.
#' @param trait_col an integer number defining the column of interested trait
#' in \code{phenofile}.
#' @param out the directory that save the results. RVFamSq produces output file
#' with the name of the gene, score, p-value of the gene, number of sample size
#' and number of families
#' @param kin a kinship matrix calculated from the pedigree using the
#' “kinship2” package
#' @param start-par an optional argument that defines the starting values of
#' the parameters fitting. User can define specific values to start the
#' parameters fitting, otherwise, the starting values will generate randomly.
#' @return \item{results}{a data frame results containing the name of the gene,
#' the value of the score and p-value write into the file named by the name of
#' the gene under the directory defined by the argument \code{out}.}
#' \item{paras}{a list of parameters estimated by maximizing the likelihood.
#' The values of the paras save to the file defined by the argument
#' \code{parafile}”. If the \code{parafile} is existed, the values of the paras
#' load from the file \code{parafile}.}
#' @author Zhihui Zhang and Suzanne Leal
#' @references Chen, W.-M., and Abecasis, G.R. (2007). Family-Based Association
#' Tests for Genomewide Association Scans. Am. J. Hum. Genet. 81, 913–926.\cr
#' @examples
#' 
#' library (kinship2)
#' library (RVFamSq)
#' 
#' ##  Define variables pointing to the path to the example data
#' data_dir<-system.file("data", package = "RVFamSq")
#' maffile<-paste0(data_dir,"/test.sfs")
#' genofile<-paste0(data_dir, "test.geno")
#' phenofile<-paste0(data_dir,"/test.pheno")
#' 
#' ##  Define output directory that save the results
#' out<-"/directory-to-save-the-results"
#' 
#' ## Define the files that save the parameters estimated under the null model.
#' ## If the file is not existed, the package will estimate the parameters based on the available data.
#' parafile<-"/ directory-to-example-data /paras.rds"
#' 
#' ## Load phenotype data and calculate the kinship matrix of the pedigree.
#' ped_pheno<-read.table(phenofile)
#' kin_pre<-data.frame(id=ped_pheno[,2], mom=ped_pheno[,4], dad=ped_pheno[,3], sex=ped_pheno[,5])
#' tped <- with(kin_pre, pedigree(id, dad, mom, sex, famid=ped_pheno[,1]))
#' kin <- kinship(tped)
#' 
#' ## Run RVFamsSq package and calculate the statistical score of the interested gene.
#' RV_FamSq(pedfile, phenofile, maffile, parafile, c(5,7), 6, out,kin)
#'
#' @import bbmle 
#' @import mvtnorm 
#' @import rlist
#' @export
RV_FamSq <- function(genofile, phenofile, maffile, parafile,covar_col, trait_col, out,kin, start_par=NULL) {
  library(bbmle)
  library(mvtnorm)
  library(rlist)

  if(!file.exists(genofile)) stop("Genotype data does not exist.")
  if(!file.exists(maffile)) stop("MAF data does not exist.")
  if(!file.exists(phenofile)) stop("Phenotype data does not exist.")


  ped_pheno<-as.matrix(read.table(phenofile))
  ped_geno<-as.matrix(read.table(genofile))
  maf_data<-as.matrix(read.table(maffile))

  if (ped_pheno[,2] %contain% ped_geno[,2]) {

    miss_index<-which(is.na(ped_pheno[,covar_col]) | is.na(ped_pheno[,trait_col]))

    if(length(miss_index)>0) {
      print("Delete samples without the data of phenotype and covariant")
      ped_pheno<-ped_pheno[-miss_index,]
    }

    ped_pheno<-ped_pheno[,c(1:4,covar_col,trait_col)]
    geno_index<-match(ped_pheno[,2],ped_geno[,2])
    ped_data<-cbind(ped_pheno,ped_geno[geno_index,3:ncol(ped_geno)])

  } else {
    stop("Some samples in phenotype file miss the genotype data")
  }

  covar_update_col<-c(5:(4+length(covar_col)))
  trait_update_col<-c(5+length(covar_col))
  geno_update_col<-c((6+length(covar_col)):ncol(ped_data))

  genos<-ped_data[,geno_update_col]

  if(missing(kin)) {
    stop("Kinship matrix is missing.")
  } else {
    if (rownames(kin) %contain% ped_data[,2]) {
      assign("kin", kin, envir=.GlobalEnv)
    } else {
      stop("Some samples in pedfile are not included in Kinship matrix")
    }
  }

  if (ncol(genos)!=nrow(maf_data)) stop("Dimension of genotype and number of varints in MAF file does not match")

  if(!file.exists(parafile)) {
    par_num<-3+length(covar_col)
    if(is.null(start_par)) {
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
    assign("simdata", list(id=ped_data[,2],pheno=ped_data[,trait_update_col],covar=ped_data[,covar_update_col],famids=ped_data[,1],genos=genos, maf=maf_data[,4]), envir=.GlobalEnv)
    print(paste("Estimating the values of parameters",paste(names(paras),collapse = ","),sep = " "))
    nullmod<-nlf_fit(paras)
    paras<-as.list(coef(nullmod))
    list.save(paras, parafile)
    RV_FamSq(pedfile, phenofile, maffile, parafile,covar_col, trait_col, out,kin)
  } else {
    print("Reading parameters from parafile... ")
    paras<-list.load(parafile)
    print(unlist(paras))

    peds_grp_by_fam<-split(seq_len(nrow(ped_data)),ped_data[,1])
    scores<-unlist(lapply(peds_grp_by_fam, function(x) Tscore(ped_data[unlist(x),],as.numeric(maf_data[,4]),paras,covar_update_col,trait_update_col,geno_update_col)))
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
