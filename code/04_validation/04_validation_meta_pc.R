#!/usr/bin/env Rscript
# Stage 4 refinement: meta-analysis for the protein-coding signature.
## --- portable project roots: do not edit (see config/paths.R) ---------------
.start_dir <- tryCatch({
  .a <- commandArgs(trailingOnly = FALSE)
  dirname(normalizePath(sub("^--file=", "", .a[grep("^--file=", .a)][1])))
}, error = function(e) getwd(), warning = function(w) getwd())
.find_root <- function(d) {
  repeat {
    if (file.exists(file.path(d, "config", "paths.R"))) return(d)
    .p <- dirname(d)
    if (identical(.p, d)) stop("repository root not found -- set GLIOMA_ROOT")
    d <- .p
  }
}
REPO_ROOT <- if (nzchar(Sys.getenv("GLIOMA_ROOT"))) normalizePath(Sys.getenv("GLIOMA_ROOT")) else .find_root(.start_dir)
source(file.path(REPO_ROOT, "config", "paths.R"))
## ----------------------------------------------------------------------------

suppressMessages({library(survival); library(metafor)})
set.seed(20260909)
setwd(ROOT)
RES <- "results"; D <- file.path(DATA_ROOT, "CGGA")
sig <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors=FALSE)
sm <- setNames(sig$coefficient, sig$gene)
cind <- function(risk,t,e){n<-length(risk);c<-0;d<-0;for(i in 1:n){if(!e[i])next;ti<-t[i]
  for(j in 1:n){if(i==j)next;if(t[j]>ti||(t[j]==ti&&!e[j])){d<-d+1;if(risk[i]>risk[j])c<-c+1 else if(risk[i]==risk[j])c<-c+0.5}}}
  if(d==0)NA else c/d}
boot <- function(risk,t,e,B=300){ci<-numeric(B);n<-length(risk)
  for(b in 1:B){i<-sample(n,replace=TRUE);ci[b]<-cind(risk[i],t[i],e[i])}
  c(mean=mean(ci,na.rm=TRUE),lo=unname(quantile(ci,.025,na.rm=TRUE)),hi=unname(quantile(ci,.975,na.rm=TRUE)),se=sd(ci,na.rm=TRUE))}

rows <- list()
for (name in c("693","325")) {
  zf <- paste0(D,"/CGGA.mRNAseq_",name,".RSEM-genes.20200506.txt.zip"); tmp <- tempfile(fileext=".txt")
  system(paste0("unzip -p '",zf,"' 'CGGA.mRNAseq_",name,".RSEM-genes.20200506.txt' > '",tmp,"'"))
  ex <- read.delim(tmp,sep="\t",check.names=FALSE); rownames(ex)<-ex[,1]; ex<-as.matrix(ex[,-1])
  cl <- read.delim(paste0(D,"/CGGA.mRNAseq_",name,"_clinical_clean.tsv"),stringsAsFactors=FALSE)
  prim<-cl[cl$PRS_type=="Primary",]; prim$OS<-as.numeric(prim$OS); prim$Censor<-as.numeric(prim$Censor)
  clp<-prim[!is.na(prim$OS)&prim$OS>0,]; ss<-intersect(colnames(ex),clp$CGGA_ID); ex<-ex[,ss]; clp<-clp[match(ss,clp$CGGA_ID),]
  X<-log2(as.matrix(ex)+1); g<-intersect(names(sm),rownames(X)); risk<-as.numeric(t(X[g,,drop=FALSE])%*%sm[g])
  bc<-boot(risk,clp$OS,clp$Censor)
  rows[[name]]<-data.frame(cohort=paste0("CGGA-",name),n=length(risk),events=sum(clp$Censor),genes=length(g),
                           C_index=bc[["mean"]],lo=bc[["lo"]],hi=bc[["hi"]],se=bc[["se"]])
  cat(sprintf("CGGA-%s: C=%.3f [%.3f-%.3f] (genes %d/%d)\n",name,bc[["mean"]],bc[["lo"]],bc[["hi"]],length(g),length(sm)))
}
df<-do.call(rbind,rows)
rf<-rma(yi=df$C_index,sei=df$se,method="REML")
df<-rbind(df,data.frame(cohort="POOLED (REML)",n=sum(df$n),events=sum(df$events),genes=NA,
                        C_index=rf$b[1],lo=rf$ci.lb,hi=rf$ci.ub,se=rf$se))
cat(sprintf("\nPOOLED protein-coding C-index = %.3f [%.3f-%.3f], I2=%.1f\n",rf$b[1],rf$ci.lb,rf$ci.ub,rf$I2))
write.table(df,file.path(RES,"external_validation_summary_pc.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("DONE\n")
