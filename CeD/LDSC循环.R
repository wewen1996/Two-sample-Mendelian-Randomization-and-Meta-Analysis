#setwd(dir="D:/mendelian randomization/test/test")  #设置你自己的工作路径

# -加载需要用到的包，没安装的需要提前安装一下
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(dplyr)


# 暴露，需要用的文件放在你设置的工作路径下
FileNames <-list.files(paste0(getwd()),pattern=".gz")
exp_dat_ids <- FileNames
exps <- FileNames


# 结局
##从网上下载vcf数据，读入后使用
library(VariantAnnotation)
data0=VariantAnnotation::readVcf("ebi-a-GCST000612.vcf.gz")
library(gwasvcf)
data1=gwasvcf::vcf_to_tibble(data0)
###缺少gwasvcf, genetics.binaRies, ldscr包，一直装不上
# out<-fread("BMI.txt",header = T)#读取本地结局（这里将文件名换成你自己本地结局的文件名）
# out$trait <- 'BMI'  #改成你所用到的结局的表型的名称
# outcomeid <- out
# rm(out)
# head(outcomeid)
# outcomeid <- c("ebi-a-GCST000612")
# exp1 <- extract_instruments(outcomes = "ebi-a-GCST000612",p1 = 1e-01,)
# exp1$trait <- 'CeD'
# out_data <- extract_outcome_data(snps= exposure$SNP,outcomes = outcomeid)
data1$trait <- 'CeD'
outcomeid <- data1
##对数据进行格式转换，注意列名
outcome<-format_data(outcomeid,type="outcome",
                     snp_col = "SNP",
                     phenotype_col = "trait",
                     effect_allele_col = "ALT",
                     other_allele_col = "REF",
                     beta_col = "beta",
                     se_col = "SE",
                     samplesize_col = "N",
                     pval_col = "p",
                     eaf_col = "AF",
                     chr_col = "chr",
                     pos_col = "pos"
)
rm(outcomeid)

#
dir.create(path = "mendelian")  #在你的工作路径下创建一个名为"mendelian"的文件夹，后面的结果都会生成在这个文件夹里面

#######以下为循环代码，不需要进行更改######################
qaq = 1
for (qaq in 1:length(exp_dat_ids)) { # 
  exp_dat_id <- exp_dat_ids[qaq]
  exp <- exps[qaq]
  
  ##逐个读入暴露数据
  exposure<- try(fread(paste0(getwd(),"/",FileNames[qaq]),fill=TRUE),silent = T)
  
  exposure$Phenotype <- FileNames[qaq]
  
  head(exposure)

###数据格式转换，注意列名
  exposure<-format_data(exposure,type="exposure",
                   snp_col = "MarkerName",
                   phenotype_col = "Phenotype",
                   effect_allele_col = "Allele1",
                   other_allele_col = "Allele2",
                   beta_col = "Effect",
                   se_col = "StdErr",
                   samplesize_col = "TotalSampleSize",
                   pval_col = "`P-value`",
                   eaf_col = "Freq1",
                   chr_col = "Chr",
                   pos_col = "Pos")
  
  
  head(exposure)
  head(outcome)
  

  exposure$id.exposure<-FileNames[qaq]
  exposure$exposure<-FileNames[qaq]
  outcome$id.outcome<-outcome$outcome
  

  ###依赖文件放到工作路径内
  ld="G:/MR/血清代谢物课程（壹隅_）/5.LDSC代码＋视频教学/依赖文件/LDSC依赖文件/1000G" 
  wld="G:/MR/血清代谢物课程（壹隅_）/5.LDSC代码＋视频教学/依赖文件/LDSC依赖文件/1000G"
  

  LDSC_rg<-function(expo,outcome,an,sample_prev=NA,
                    population_prev=NA,ld,wld,chr_filter=c(1:22),n_blocks=200){
    id.o<-outcome$id.outcome[1]
    id.e<-expo$id.exposure[1]
    
    expo<-expo%>%mutate(Z=beta.exposure/se.exposure)
    expo<-expo%>%select(SNP=SNP,N=samplesize.exposure,Z=Z
                        ,A1=effect_allele.exposure
                        ,A2=other_allele.exposure)
    expo<-as_tibble(expo)
    
    outcome<-outcome%>%mutate(Z=beta.outcome/se.outcome)
    outcome<-outcome%>%select(SNP=SNP,N=samplesize.outcome,Z=Z
                              ,A1=effect_allele.outcome
                              ,A2=other_allele.outcome)
    outcome<-as_tibble(outcome)
    
    
    dat<-list(expo,outcome)
    names(dat)<-c(id.e,id.o)
    
    rm(expo,outcome)
    
    
    res<-try(ldscr::ldsc_rg(dat,ancestry = an,sample_prev=sample_prev,
                            population_prev=population_prev,ld=ld,wld=wld,
                            n_blocks=n_blocks,chr_filter=chr_filter))
    
    return(res)
    
  }
  
  LDSC_res<-LDSC_rg(exposure,outcome,ld=ld,wld=wld)
  
  h2 <- LDSC_res[["h2"]]
  rg <- LDSC_res[["rg"]]
  
  
  
  write.csv(h2,file = paste0("mendelian/",exp,"_h2.csv"), row.names = FALSE)
  write.csv(rg,file = paste0("mendelian/",exp,"_rg.csv"), row.names = FALSE)
  
  }

##结果合并
fs=list.files("D:/mendelian randomization/test/test/mendelian", pattern = "rg.csv",full.names = TRUE) #把导出的路径改成自己的工作路径  
df = map_dfr(fs, read.csv)
write.csv(df,"rg_res.csv")

fs=list.files("D:/mendelian randomization/test/test/mendelian", pattern = "h2.csv",full.names = TRUE) #把导出的路径改成自己的工作路径  
df = map_dfr(fs, read.csv)
write.csv(df,"h2_res.csv")
                        