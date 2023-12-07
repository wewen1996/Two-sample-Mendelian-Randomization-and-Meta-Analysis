#setwd(dir="D:/mendelian randomization/test")  #设置你自己的工作路径

# -加载需要用到的包，没安装的需要提前安装一下
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(ieugwasr)

# 暴露，把文件放在你设置的工作路径下
FileNames <-list.files(paste0(getwd()),pattern=".gz")
exp_dat_ids <- FileNames
exps <- FileNames

# 结局,把编号改成自己需要的结局的编号
# 结局
out<-fread("IA_stage_1.txt",header = T)#读取本地结局（这里将文件名换成你自己本地结局的文件名）
out$trait <- 'IA'  #改成你所用到的结局的表型的名称
outcomeid <- out
rm(out)
head(outcomeid)
#
dir.create(path = "mendelian 88")  #在你的工作路径下创建一个名为"mendelian 88 IA-YZ"的文件夹，后面的结果都会生成在这个文件夹里面

#######以下为循环代码，不需要进行更改######################
qaq <- 1
#for (qaq in 1:length(exp_dat_ids)) { # 
  exp_dat_id <- exp_dat_ids[qaq]
  exp <- exps[qaq]
  
  d3<- try(fread(paste0(getwd(),"/",FileNames[qaq]),fill=TRUE),silent = T)
  d3<-subset(d3,d3$`P-value`<1e-5)
  # rm(d1)
  # d3<-d2[,c(1,2,3,4,8,9,10,16)]
  d3$PHENO <- FileNames[qaq]
  names(d3)[names(d3) == 'MarkerName'] <- 'SNP'
  
  
 d3<-format_data(d3,
                 type="exposure",
                 phenotype_col = "PHENO",
                 snp_col = "SNP",
                 beta_col = "Effect",
                 se_col = "StdErr",
                 pval_col = "P-value",
                 samplesize_col = "TotalSampleSize",
                 eaf_col = "Freq1",
                 effect_allele_col = "Allele1",
                 other_allele_col = "Allele2")
 
#exp_data <- clump_data(d3,clump_kb = 500,clump_r2 = 0.01)
     
 library(ieugwasr)     
  
 d4<- ld_clump(
    #dat = X1,
    clump_kb = 500,
    clump_r2 = 0.01,
    pop = "EUR",
    dplyr::tibble(rsid=d3$SNP, pval=d3$pval.exposure, id=d3$id.exposure),
    #get_plink_exe()
    plink_bin = "E:/R-4.3.0/library/plinkbinr/bin/plink_Windows.exe",
    #欧洲人群参考基因组位置
    bfile = "D:/mendelian randomization/test1/依赖文件/本地LD依赖文件及MAF文件/1kg.v3/1kg.v3/EUR"
  )
  
  exp_data<-subset(d3,SNP %in% d4$rsid) 
  
  if(length(exp_data[,1])>0){
    outcome_dat<-merge(exp_data,outcomeid,by.x = "SNP",by.y = "SNP")#修改本地结局数据中对应rsid列的列名为SNP
    write.csv(outcome_dat,file = "d.csv")
    out_data <- read_outcome_data(
      snps = exp_data$SNP, 
      filename = "d.csv",
      sep = ",",
      phenotype_col = "trait",
      snp_col = "SNP",
      beta_col = "be",
      se_col = "se",
      eaf_col="freq",
      samplesize_col = "n",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      pval_col = "p")
    
    
    if(length(out_data[,1])>0){  
      dat <- TwoSampleMR::harmonise_data(
        exposure_dat = exp_data,
        outcome_dat = out_data)
      
      ####回文的直接去除
      dat <-subset(dat,mr_keep==TRUE)
      
     dat$samplesize.outcome <- 300000
       
      #计算F值和R2
      get_f<-function(dat,F_value=10){
        log<-is.na(dat$eaf.exposure)
        log<-unique(log)
        if(length(log)==1)
        {if(log==TRUE){
          print("数据不包含eaf，无法计算F统计量")
          return(dat)}
        }
        if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("数据不包含beta，无法计算F统计量")
          return(dat)}
        if(is.null(dat$se.exposure[1])==T || is.na(dat$se.exposure[1])==T){print("数据不包含se，无法计算F统计量")
          return(dat)}
        if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("数据不包含samplesize(样本量)，无法计算F统计量")
          return(dat)}
        
        
        if("FALSE"%in%log && is.null(dat$beta.exposure[1])==F && is.na(dat$beta.exposure[1])==F && is.null(dat$se.exposure[1])==F && is.na(dat$se.exposure[1])==F && is.null(dat$samplesize.exposure[1])==F && is.na(dat$samplesize.exposure[1])==F){
          R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
          F<- (dat$samplesize.exposure-2)*R2/(1-R2)
          dat$R2<-R2
          dat$F<-F
          dat<-subset(dat,F>F_value)
          return(dat)
        }
      }
      
      dat <- get_f(dat, F_value = 10)
      
      
      res=TwoSampleMR::mr(dat,method_list= c("mr_ivw" ,
                                             "mr_weighted_median" ,
                                             "mr_egger_regression",
                                             "mr_simple_mode",
                                             "mr_weighted_mode",
                                             "mr_wald_ratio"))
      
      print(paste0(exp,"_SNP数_",res$nsnp[1]))
      
      results <- TwoSampleMR::generate_odds_ratios(res)
      
      results$estimate <- paste0(
        format(round(results$or, 2), nsmall = 2), " (", 
        format(round(results$or_lci95, 2), nsmall = 2), "-",
        format(round(results$or_uci95, 2), nsmall = 2), ")")
      resdata <- dat
      openxlsx::write.xlsx(dat,file = paste0("mendelian 88/",exp,"-dat.xlsx"), row.names = FALSE)
      
      names(resdata)
      Assumption13 <- subset(resdata,mr_keep==TRUE,
                             select = c("SNP","pval.exposure",
                                        "pval.outcome", # "F_statistic",
                                        "mr_keep"))
      
      openxlsx::write.xlsx(x = list(
        "main"=results,
        "Assumption13"=Assumption13),
        overwrite = TRUE,
        paste0("mendelian 88/",exp,"-res.xlsx"))
      
    }}
  if(length(dat[,1])>2){
    res_hete <- TwoSampleMR::mr_heterogeneity(dat)
    res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
    res_leaveone <- mr_leaveoneout(dat)  # 
    
 
         #若需运行MR—presso,将最左侧7个#删去即可 
    #         res_presso <- TwoSampleMR::run_mr_presso(dat,
    #                                             NbDistribution = 100)
           # [["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
    #       sink(paste0("mendelian 88/",exp,"_PRESSO.txt"),
    #                        append=FALSE,split = FALSE) 
    #       print(res_presso)
    #       sink()
    #       print(res_presso)
    
    
    
    
    p1 <- mr_scatter_plot(res, dat)
    p1[[1]]
    pdf(paste0("mendelian 88/",exp,"_scatter.pdf"))
    print(p1[[1]])
    dev.off()
    
    res_single <- mr_singlesnp(dat)
    p2 <- mr_forest_plot(res_single)
    pdf(paste0("mendelian 88/",exp,"_forest.pdf"))
    print(p2[[1]])
    dev.off()
    
    p3 <- mr_funnel_plot(res_single)
    pdf(paste0("mendelian 88/",exp,"_funnel.pdf"))
    print(p3[[1]])
    dev.off()
    
    res_loo <- mr_leaveoneout(dat)
    pdf(paste0("mendelian 88/",exp,"_leave_one_out.pdf"))
    print(mr_leaveoneout_plot(res_loo))
    dev.off()
    
    
    library(magrittr)
    res3 <- res[1:3,]
    
    
    
    # 转换成论文格式
    library(magrittr)
    # Main result 
    res4 <- tidyr::pivot_wider(
      res3,names_from ="method",names_vary = "slowest",
      values_from = c("b","se","pval") )
    # Heterogeneity statistics
    res_hete2 <- tidyr::pivot_wider(
      res_hete,names_from ="method",names_vary = "slowest",
      values_from = c("Q","Q_df","Q_pval") ) %>% 
      dplyr::select( -id.exposure,-id.outcome,-outcome,-exposure)
    # Horizontal pleiotropy
    res_plei2 <- dplyr::select(res_plei,
                               egger_intercept,se,pval)
  
    
    # Merge
    res_ALL <- cbind(res4, res_hete2, res_plei2)
    
    write.csv(res_ALL,file = paste0("mendelian 88/",exp,".csv"), row.names = FALSE)
  }}


#把导出的路径改成自己的工作路径  
fs=list.files("工作路径/mendelian 88", pattern = "csv",full.names = TRUE) #把导出的路径改成自己的工作路径  
df = map_dfr(fs, read.csv)
write.csv(df,"res.csv")
