setwd("./")

library(data.table)
library(dplyr)
library(ggplot2)
library(broom)
library(pROC)
library(gtsummary)
library(MatchIt)
library(survival)
library(casebase)

TPMI_pheno_table<-fread("./pheno/TPMI_phencode_pheno.processed.txt",sep="\t")
TPMI_imp<-fread("./data/PRSCatlog-ALL-done-imputed-230310.txt",sep="\t")

### Add PC 230901 ####

PGS_phecode<-fread("./phencode/phecode_PGS_codebook_addPC_0901.txt",sep="\t")
PGS_phecode<-PGS_phecode[,c("phenotype","PGSID")]

### PC table ####

PC_table<-fread("./PC/CMUH_30W_PCs.txt",sep="\t")

################

ready_process_list <- fread("./pheno/Ready_list_done_addPC_230901.txt",sep="\t")
ready_process_list_ready<-ready_process_list[c(1:114),]
ready_process_list_done<-ready_process_list_ready$PhenotypeName

ratios<-c("output-1to2_matching","output-1to4_matching","output-1to6_matching","output-1to8_matching")

# Ready to process

dir.create(paste0("./output/"))

for (t in 1:length(ratios)){

  dir.create(paste0("./output/",ratios[t],"/"))

  for (i in 1:length(ready_process_list_done)){
    
 
    interest_list<-PGS_phecode[grepl(ready_process_list_done[i],PGS_phecode$phenotype),]
    interest_list<-subset(PGS_phecode,PGS_phecode$phenotype==ready_process_list_done[i])
    interest_PGS<-interest_list$PGSID
    
    Full_PGS_list<-c("IID",interest_PGS)
    PGS_interested<-select(TPMI_imp, matches(Full_PGS_list))
    
    pattern = "PGS"
    pattern_PGS <- PGS_interested[,grep(pattern = pattern, colnames(PGS_interested))]
    list_PGS<-colnames(PGS_interested)[pattern_PGS]
    
    PRS.result.df<-left_join(TPMI_pheno_table,PGS_interested,by=c("GenotypeName"="IID"))
    PGS_NUM <- list_PGS 

    dir.create(paste0("./output/",ratios[t],"/",ready_process_list_done[i]))
    
    Pheno_col<-ready_process_list_done[i]     
        
    for (j in 1:length(PGS_NUM)){
            
       folder <- dir.create(paste0("./output/",ratios[t],"/",Pheno_col,"/",PGS_NUM[j]))
       sub_folder <- dir.create(paste0("./output/",ratios[t],"/",Pheno_col,"/",PGS_NUM[j],'/','plot'))

    ########################################################
    cat("\n")
    cat("###############################################")
    cat("\n")
    cat("\n")
    cat("The Program is now processing:",ready_process_list_done[i])
    cat("\n")
    cat("\n")
    cat("###############################################")
    cat("\n")
    cat("\n")
    cat("The Program is now processing:",PGS_NUM[j])
    cat("\n")
    cat("\n")
    cat("###############################################\n")
  #########################################################
  
        
 PRS.result.plot.tmp <- paste0("./Result/Other_matching_ratio_230906/Ratios/",ratios[t],"/",ready_process_list_done[i],"/",PGS_NUM[j],"/plot/",PGS_NUM[j],".plot-data.txt",collapse = '')
 
 PRS.result.plot<-fread(PRS.result.plot.tmp,sep="\t")

 PRS.result.plot$Pheno<-recode_factor(PRS.result.plot$Pheno,"Control"="0","Case"="1")

######### Prepare model data ########

  PRS.result.model<-PRS.result.plot

############################### Ready GLM PRS Model Data (Drop NA's) ########################################

  PRS.pheno.model=PRS.result.model
  PRS.pheno.model=na.omit(PRS.pheno.model)

  colnames(PRS.pheno.model)<-c("PatientID","IID","Sex","Age","Pheno","SCORE","PC1","PC2","PC3","PC4","percentile")
  
  # To avoid only one sex error
  PRS.pheno.model$Sex<-recode_factor(PRS.pheno.model$Sex,"Male"="2","Female"="1")
  PRS.pheno.model$Sex<-as.numeric(PRS.pheno.model$Sex)

  cat('\n')
  cat('\n')
  cat('\n ###########################################')
  cat('\n')
  cat('\n')
  cat('\n PRS GLM model has dropped row with NAs....')
  cat('\n')
  cat('\n')
  cat('\n ###########################################')

##### 80% of the sample size #### 

  smp_size <- floor(0.8 * nrow(PRS.pheno.model))

## set the seed to make your partition reproducible
  set.seed(123)
  train_ind <- sample(seq_len(nrow(PRS.pheno.model)), size = smp_size)

##### Train data #####
  PRS.train.pheno <- PRS.pheno.model[train_ind, ]

##### Test data #####
  PRS.test.pheno <- PRS.pheno.model[-train_ind, ]

##### Train Test data output #####

  tmp.tarin.data<-paste0("./output/",ratios[t],"/",Pheno_col,"/",PGS_NUM[j],'/plot/',PGS_NUM[j],'.AUC.train-data.txt',collapse = '')
  fwrite(PRS.train.pheno,tmp.tarin.data,sep="\t",col.names=T)

  tmp.test.data<-paste0("./output/",ratios[t],"/",Pheno_col,"/",PGS_NUM[j],'/plot/',PGS_NUM[j],'.AUC.test-data.txt',collapse = '')
  fwrite(PRS.test.pheno,tmp.test.data,sep="\t",col.names=T)

### Summary of data ####
  cat('\n')
  cat('All data of PRS caculation:')
  cat('\n')
  cat('\n')
  summary(PRS.pheno.model)
  cat('\n')
  cat('\n')
  cat('Train data of PRS caculation:')
  cat('\n')
  cat('\n')
  summary(PRS.train.pheno)
  cat('\n')
  cat('\n')
  cat('Test data of PRS caculation:')
  cat('\n')
  cat('\n')
  summary(PRS.test.pheno)
  cat('\n')
  cat('\n')

################################ Read GLM PRS Model Data ###############################
  
  mod_0 <- glm( Pheno ~ Age+Sex, data=PRS.train.pheno, family="binomial")
  mod_1 <- glm( Pheno ~ SCORE, data=PRS.train.pheno, family="binomial")
  mod_2 <- glm( Pheno ~ SCORE+Age+Sex, data=PRS.train.pheno, family="binomial")
  mod_3 <- glm( Pheno ~ SCORE+Age+Sex+PC1+PC2+PC3+PC4, data=PRS.train.pheno, family="binomial")

##################################### Model 0 (Sex Age only) #####################################

  library(pROC)

  ## Train
  train_0_prob = predict(mod_0, data=PRS.train.pheno ,type='response')
  train_0_roc = roc(PRS.train.pheno$Pheno ~ train_0_prob, plot = FALSE, print.auc = TRUE, legacy.axes=TRUE)
  
  ## Test
  test_0_prob = predict(mod_0, newdata = PRS.test.pheno, type = "response")
  test_0_roc = roc(PRS.test.pheno$Pheno ~ test_0_prob, plot = FALSE, print.auc = TRUE, legacy.axes=TRUE)

##################################### Model 1 (Base only) #####################################


  ## Train
  train_1_prob = predict(mod_1, data=PRS.train.pheno ,type='response')
  train_1_roc = roc(PRS.train.pheno$Pheno ~ train_1_prob, plot = FALSE, print.auc = TRUE, legacy.axes=TRUE)
  
  ## Test
  test_1_prob = predict(mod_1, newdata = PRS.test.pheno, type = "response")
  test_1_roc = roc(PRS.test.pheno$Pheno ~ test_1_prob, plot = FALSE, print.auc = TRUE, legacy.axes=TRUE)


##################################### Model 2 (PRS only) #####################################

  ## Train
  train_2_prob = predict(mod_2, data=PRS.train.pheno ,type='response')
  train_2_roc = roc(PRS.train.pheno$Pheno ~ train_2_prob, plot = FALSE, print.auc = TRUE, legacy.axes=TRUE)
  
  ## Test
  test_2_prob = predict(mod_2, newdata = PRS.test.pheno, type = "response")
  test_2_roc = roc(PRS.test.pheno$Pheno ~ test_2_prob, plot = FALSE, print.auc = TRUE, legacy.axes=TRUE)

##################################### Model 3 (Full model) #####################################


  ## Train
  train_3_prob = predict(mod_3, data=PRS.train.pheno ,type='response')
  train_3_roc = roc(PRS.train.pheno$Pheno ~ train_3_prob, plot = FALSE, print.auc = TRUE, legacy.axes=TRUE)
  
  ## Test
  test_3_prob = predict(mod_3, newdata = PRS.test.pheno, type = "response")
  test_3_roc = roc(PRS.test.pheno$Pheno ~ test_3_prob, plot = FALSE, print.auc = TRUE, legacy.axes=TRUE)

######################## AUC plot Table ################################

## Table Output

#### Table 0 ######

  out_a_0=coords(test_0_roc, "best", ret = c("auc","threshold", "specificity", "sensitivity", "accuracy","precision", "recall"), transpose = FALSE, print.auc = TRUE)

  #Get the first row to prevent error of 0 rows of metrics
  out_a_0=out_a_0[1,]
  out_b_0=as.data.frame(auc(test_0_roc))

  colnames(out_b_0)<-"AUC"
  out_final_0<-cbind(out_a_0,out_b_0)

  out_final_0$Model<-c("Age + Sex model")


#### Table 1 ######

  out_a_1=coords(test_1_roc, "best", ret = c("auc","threshold", "specificity", "sensitivity", "accuracy","precision", "recall"), transpose = FALSE, print.auc = TRUE)

  #Get the first row to prevent error of 1 rows of metrics
  out_a_1=out_a_1[1,]
  out_b_1=as.data.frame(auc(test_1_roc))

  colnames(out_b_1)<-"AUC"
  out_final_1<-cbind(out_a_1,out_b_1)

  out_final_1$Model<-c("PRS model")


#### Table 2 ######

  out_a_2=coords(test_2_roc, "best", ret = c("auc","threshold", "specificity", "sensitivity", "accuracy","precision", "recall"), transpose = FALSE, print.auc = TRUE)

  #Get the first row to prevent error of 2 rows of metrics
  out_a_2=out_a_2[1,]
  out_b_2=as.data.frame(auc(test_2_roc))

  colnames(out_b_2)<-"AUC"
  out_final_2<-cbind(out_a_2,out_b_2)

  out_final_2$Model<-c("PRS + Age + Sex")

#### Table 3 ######

  out_a_3=coords(test_3_roc, "best", ret = c("auc","threshold", "specificity", "sensitivity", "accuracy","precision", "recall"), transpose = FALSE, print.auc = TRUE)

  #Get the first row to prevent error of 3 rows of metrics
  out_a_3=out_a_3[1,]
  out_b_3=as.data.frame(auc(test_3_roc))

  colnames(out_b_3)<-"AUC"
  out_final_3<-cbind(out_a_3,out_b_3)

  out_final_3$Model<-c("Full model (PRS + Age + Sex + PCs)")

########################

  out_final<-rbind(out_final_0,out_final_1,out_final_2,out_final_3)
  
  out_final$PGSID <- PGS_NUM[j]

  out_final<-out_final[,c(9,8,7,1:6)]

  out_final$PhenotypeName<-Pheno_col

  out_final<-out_final[,c(10,1:9)]

  tmp_out<-paste0("./output/",ratios[t],"/",Pheno_col,"/",PGS_NUM[j],'/plot/',PGS_NUM[j],'.Performance.txt')
  
  fwrite(out_final,tmp_out,sep="\t",col.names = T)

  ### Model plot ###

    tmp_model_plot<-paste0("./output/",ratios[t],"/",Pheno_col,"/",PGS_NUM[j],'/plot/',PGS_NUM[j],'.AUCs.png',collapse = '')

    png(tmp_model_plot,height = 500,width  = 500)

    test_0_roc = roc(PRS.test.pheno$Pheno ~ test_0_prob)
    test_1_roc = roc(PRS.test.pheno$Pheno ~ test_1_prob)
    test_2_roc = roc(PRS.test.pheno$Pheno ~ test_2_prob)
    test_3_roc = roc(PRS.test.pheno$Pheno ~ test_3_prob)

    plot(test_0_roc,print.auc = TRUE, print.auc.y = .5 , col='darkmagenta')
    plot(test_1_roc,print.auc = TRUE, print.auc.y = .4 , col='black',add=TRUE)
    plot(test_2_roc,print.auc = TRUE, print.auc.y = .3 , col='red'  ,add=TRUE)
    plot(test_3_roc,print.auc = TRUE, print.auc.y = .2 , col='blue' ,add=TRUE)
    
    text(0.15, .5, paste("Age + Sex model"))
    text(0.15, .4, paste("PRS model"))
    text(0.15, .3, paste("PRS + Age + Sex"))
    text(0.15, .2, paste("PRS + Age + Sex + PCs"))

    cat('\n')
    cat('\nAUC of PRS Model is:', auc(test_1_roc))
    cat('\n')
    cat('\nAUC of PRS + Age + Sex Model is:', auc(test_2_roc))
    cat('\n')
    cat('\nAUC of Full Model is:', auc(test_3_roc))
    dev.off()

    }
    
  }

}