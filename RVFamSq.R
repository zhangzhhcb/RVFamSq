"%contain%" <- function(values,x) {
  tx <- table(x)
  tv <- table(values)
  z <- tv[names(tx)] - tx
  all(z >= 0 & !is.na(z))
}

#define negative likelihood function
nlf <- function(miu,betax1,sitag,sitae) {
  sum_likelihood=0
  for (n in unique(simdata$famid)) {
    sid=which(simdata$famids==n)
    sid_count=length(sid)
    EYi=miu+betax1*simdata$covar[sid]
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
  m1 <- mle2(nlf,start=start, data=simdata, method = "Nelder-Mead", skip.hessian = FALSE, control=list(trace=TRUE, maxit=500))
  m2 <- mle2(nlf,start=as.list(coef(m1)), control=list(trace=TRUE, maxit=500, parscale=coef(m1)), data=simdata)
  return(m2)
}

Tscore<-function(ped,maf,paras) {
  miu<-paras$miu
  betax1<-paras$betax1
  sitag<-paras$sitag
  sitae<-paras$sitae

  EYi<-miu+betax1*ped[,5]
  yi<-ped[,6]
  kindex<- match(ped[,2],rownames(kin))

  fanij<-as.matrix(2*kin[kindex,kindex]*sitag^2)
  diag(fanij)<-sitag^2 + sitae^2
  fanij_inv<-solve(fanij)

  g<-rowSums(ped[,7:ncol(ped)]) - sum(2*maf)

  upperdot = t(g) %*% fanij_inv %*% (yi - EYi)
  lowerdot = t(g) %*% fanij_inv %*% g

  return(c(upperdot,lowerdot))
}

RV_FamSq <- function(pedfile, maffile, parafile,out,kin, start_par=NULL) {
  library(bbmle)
  library(mvtnorm)
  library(rlist)

  if(!file.exists(pedfile)) stop("Pedfile does not exist.")
  if(!file.exists(maffile)) stop("MAF-file does not exist.")

  ped_geno<-as.matrix(read.table(pedfile))
  maf_data<-as.matrix(read.table(maffile))

  miss_index<-which(is.na(ped_geno[,6]) | is.na(ped_geno[,5]))

  if(length(miss_index)>0) {
    print("Delete samples without the data of phenotype and covariant")
    ped_geno<-ped_geno[-miss_index,]
  }

  nl<-ncol(ped_geno)
  nr<-nrow(ped_geno)
  genos<-ped_geno[,7:nl]

  if(missing(kin)) {
    stop("Kinship matrix is missing.")
  } else {
    if (rownames(kin) %contain% ped_geno[,2]) {
      assign("kin", kin, envir=.GlobalEnv)
    } else {
      stop("Some samples in pedfile are not included in Kinship matrix")
    }
  }

  if (ncol(genos)!=nrow(maf_data)) stop("Dimensions of genodata and MAF-data do not match")

  if(!file.exists(parafile)) {
    if(is.null(start_par)) {
      print("Starting values of parameters are not defined, generating random initial values now...")
      start_par<-runif(4,-1,1)
    } else {
      if (length(start_par)<4) {
        print("Defined less parameters than required by the model, generating random initial values now...")
        start_par<-runif(4,-1,1)
      }
    }

    paras <- list(miu=start_par[1], betax1=start_par[2], sitag=start_par[3], sitae=start_par[4])
    assign("simdata", list(id=ped_geno[,2],pheno=ped_geno[,6],covar=ped_geno[,5],famids=ped_geno[,1],genos=genos, maf=maf_data[,4]), envir=.GlobalEnv)
    nullmod<-nlf_fit(paras)
    paras<-as.list(coef(nullmod))
    list.save(paras, parafile)
  } else {
    print("Reading parameters from parafile... ")
    paras<-list.load(parafile)
    miu<-paras$miu
    betax1<-paras$betax1
    sitag<-paras$sitag
    sitae<-paras$sitae

    peds<-split(seq_len(nrow(ped_geno)),ped_geno[,1])
    scores<-unlist(lapply(peds, function(x) Tscore(ped_geno[unlist(x),],as.numeric(maf_data[,4]),paras)))
    upper<-sum(scores[seq(1,length(scores),2)])
    lower<-sum(scores[seq(2,length(scores),2)])
    score<-upper^2/lower
    p=pchisq(score,df=1)

    gene_name=as.character(maf_data[1,1])
    out<-paste(c(out,"/", gene_name,".out"),collapse ="")
    pout<-data.frame(gene=gene_name, score=score, p=1-p)
    write.table(pout,out, row.names = FALSE, col.names = FALSE,quote=FALSE,sep="\t")
    print(pout)
  }
  return()
}
