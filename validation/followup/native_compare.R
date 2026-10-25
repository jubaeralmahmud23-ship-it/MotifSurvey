# Our graph adapter around unmodified installed sampling/survey APIs.
# JSON I/O and package/process loading are outside timing. Enumeration included.
args <- commandArgs(trailingOnly=TRUE)
.libPaths(c(args[1],.libPaths()))
suppressPackageStartupMessages(library(sampling))
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(Matrix))
library(jsonlite)

motifs <- function(d) {
  n <- d$n; ed <- d$edges+1L
  adj <- lapply(seq_len(n),function(i)integer(0))
  rows <- list(); types <- integer(0)
  for (i in seq_len(nrow(ed))) {
    u <- ed[i,1]; v <- ed[i,2]
    adj[[u]]<-c(adj[[u]],v);adj[[v]]<-c(adj[[v]],u)
    rows[[length(rows)+1L]]<-sort(c(u,v));types<-c(types,1L)
  }
  for (u in seq_len(n)) {
    a<-adj[[u]]
    if(length(a)>1) {
      pairs<-combn(a,2)
      for(j in seq_len(ncol(pairs))) {
        rows[[length(rows)+1L]]<-sort(c(u,pairs[,j]));types<-c(types,2L)
      }
    }
    for(v in a[a>u]) {
      for(w in intersect(a,adj[[v]])) if(w>v) {
        rows[[length(rows)+1L]]<-c(u,v,w);types<-c(types,3L)
      }
    }
  }
  m<-length(rows)
  incidence<-sparseMatrix(i=rep(seq_len(m),lengths(rows)),j=unlist(rows),x=1,dims=c(m,n))
  Y<-matrix(0,m,3);Y[cbind(seq_len(m),types)]<-1
  list(A=incidence,Y=Y,pi=vapply(rows,function(s)prod(d$p[s]),numeric(1)))
}

run_one <- function(d,method) {
  start<-proc.time()[3]
  z<-motifs(d);m<-nrow(z$A)
  after_enum<-proc.time()[3]
  # overlap product is 1 for disjoint motif units, whose weighted covariance is zero.
  logshared<-z$A %*% Diagonal(x=log(d$p)) %*% t(z$A)
  if(method=='sampling') {
    joint<-outer(z$pi,z$pi)/exp(as.matrix(logshared))
    after_setup<-proc.time()[3]
    V<-matrix(0,3,3)
    for(j in 1:3)V[j,j]<-sampling::varHT(z$Y[,j],joint,method=1)
    for(j in 1:2)for(k in (j+1):3) {
      V[j,k]<-(sampling::varHT(z$Y[,j]+z$Y[,k],joint,method=1)-V[j,j]-V[k,k])/2
      V[k,j]<-V[j,k]
    }
    totals<-colSums(z$Y/z$pi)
  } else {
    logshared@x <- -expm1(logshared@x)
    data<-data.frame(a=z$Y[,1],b=z$Y[,2],c=z$Y[,3],pi=z$pi)
    des<-survey::svydesign(ids=~1,probs=~pi,fpc=~pi,data=data,pps=survey::ppscov(logshared,weighted=TRUE),variance='HT')
    after_setup<-proc.time()[3]
    result<-survey::svytotal(~a+b+c,des)
    totals<-as.numeric(coef(result));V<-unname(vcov(result))
  }
  end<-proc.time()[3]
  list(name=d$name,method=method,units=m,totals=totals,covariance=V,enumeration_seconds=unname(after_enum-start),setup_seconds=unname(after_setup-after_enum),inference_seconds=unname(end-after_setup),total_seconds=unname(end-start))
}
input<-jsonlite::fromJSON(args[2],simplifyVector=FALSE)
results<-list()
for(d in input) {
  d$n<-as.integer(d$n);d$p<-unlist(d$p);d$edges<-matrix(unlist(d$edges),ncol=2,byrow=TRUE)
  for(method in d$methods)for(repeat_id in seq_len(d$repeats)) {
    r<-run_one(d,method);r$repeat_id<-repeat_id;results[[length(results)+1L]]<-r
    cat(d$name,method,repeat_id,r$total_seconds,'seconds\n');flush.console()
  }
}
jsonlite::write_json(list(R=R.version.string,versions=as.list(sapply(c('sampling','survey','Matrix','jsonlite'),function(p)as.character(packageVersion(p)))),results=results),args[3],auto_unbox=TRUE,digits=17,pretty=TRUE)
