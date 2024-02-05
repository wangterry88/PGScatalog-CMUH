setwd("/media/volume1/Terry/Project/CMUH-PRS-Phecode")

##### transport ######
if("survival" %in% rownames(installed.packages()) == FALSE) {
        install.packages("survival",repos = "http://cran.us.r-project.org")
}
if("casebase" %in% rownames(installed.packages()) == FALSE) {
        install.packages("casebase",repos = "http://cran.us.r-project.org")
}

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
TPMI_imp<-fread("./data/TPMI-imputed/PRSCatlog-ALL-done-imputed-230310.txt",sep="\t")

### Add PC 230901 ####

PGS_phecode<-fread("./phencode/Phencode_PGS_codebook/phecode_PGS_codebook_addPC_0901.txt",sep="\t")
PGS_phecode<-PGS_phecode[,c("phenotype","PGSID")]

### PC table ####

PC_table<-fread("./PC/CMUH_30W_PCs.txt",sep="\t")

################

ready_process_list <- fread("./pheno/Ready_list_done_addPC_230901.txt",sep="\t")
ready_process_list_ready<-ready_process_list[c(158:180),]
ready_process_list_done<-ready_process_list_ready$PhenotypeName

# Ready to process

dir.create(paste0("./output/"))

for (i in 1:length(ready_process_list_done)){
    
 
    #interest_list<-PGS_phecode[grepl(ready_process_list_done[i],PGS_phecode$phenotype),]
    interest_list<-subset(PGS_phecode,PGS_phecode$phenotype==ready_process_list_done[i])
    interest_PGS<-interest_list$PGSID
    
    Full_PGS_list<-c("IID",interest_PGS)
    PGS_interested<-select(TPMI_imp, matches(Full_PGS_list))
    
    pattern = "PGS"
    pattern_PGS <- PGS_interested[,grep(pattern = pattern, colnames(PGS_interested))]
    list_PGS<-colnames(PGS_interested)[pattern_PGS]
    
    PRS.result.df<-left_join(TPMI_pheno_table,PGS_interested,by=c("GenotypeName"="IID"))
    PGS_NUM <- list_PGS 

    dir.create(paste0("./output/",ready_process_list_done[i]))
    
    Pheno_col<-ready_process_list_done[i]     
        
    for (j in 1:length(PGS_NUM)){
            
       folder <- dir.create(paste0("./output/",Pheno_col,"/",PGS_NUM[j]))
       sub_folder <- dir.create(paste0("./output/",Pheno_col,"/",PGS_NUM[j],'/','plot'))

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
  
  select<-PGS_NUM[j]
  select_column<-c("PatientID","GenotypeName","Sex","Age",Pheno_col,select)
  
  PRS.result<-select(PRS.result.df,one_of(select_column))
  colnames(PRS.result)[5]="Pheno"

  # Add PC information to the table

  PRS.result<-left_join(PRS.result,PC_table,by=c("GenotypeName"="IID"))
  PRS.result<-na.omit(PRS.result)
  
######### Perform matching ##########
  head(PRS.result)
  PRS.result.match      <- matchit(Pheno ~ Age+Sex, data = PRS.result, method="nearest", ratio=4)
  PRS.result.match.df   <- match.data(PRS.result.match)   
  PRS.result.match.done <- PRS.result.match.df[,c(1:10)]

  #print(PRS.result.match.df)

######### Prepare model data ########  
  
  #PRS.result.plot<-PRS.result     
  PRS.result.plot<-PRS.result.match.done
  PRS.result.plot$Sex<-recode_factor(PRS.result.plot$Sex,"1"="Male","2"="Female")
  PRS.result.plot$Pheno<-recode_factor(PRS.result.plot$Pheno,"0"="Control","1"="Case")
  #print(PRS.result.plot)  
        
######### Prepare model data ########

  PRS.result.model<-PRS.result.plot

########## Ready plot data #########

  PRS.pheno.plot<-PRS.result.plot
  colnames(PRS.pheno.plot)<-c("PatientID","IID","Sex","Age","Pheno","SCORE","PC1","PC2","PC3","PC4")

# Plot data prepare
  PRS.pheno.plot$percentile<-ntile(PRS.pheno.plot$SCORE,10)

# Plot data output
  tmp_PRS.pheno.plot<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.plot-data.txt',collapse = '')
  fwrite(PRS.pheno.plot,tmp_PRS.pheno.plot,sep="\t",col.names=T)

# Plot info

  PRS.pheno.high10<-subset(PRS.pheno.plot,PRS.pheno.plot$percentile==10)
  PRS.pheno.low10<-subset(PRS.pheno.plot,PRS.pheno.plot$percentile==1)

  high10<-subset(PRS.pheno.plot,PRS.pheno.plot$percentile=="10")
  high10_line<-min(high10$SCORE)

  low10<-subset(PRS.pheno.plot,PRS.pheno.plot$percentile=="1")
  low10_line<-max(low10$SCORE)

  PRS.pheno.plot$percentile<-as.factor(PRS.pheno.plot$percentile)

  head(PRS.pheno.plot)
##################### Start Plotting ##############################

########## Distribution plot ##########

  tmp_dist_plot<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.distribution.png',collapse = '')
  

  dist_plot =ggplot(PRS.pheno.plot, aes(x=SCORE, fill=Pheno)) +
              geom_vline(aes(xintercept=high10_line), colour="#BB0000", linetype="dashed")+
              geom_vline(aes(xintercept=low10_line), colour="#BB0000", linetype="dashed")+
              geom_text(aes(x = high10_line,y=0.4,label = "High 10% PRS"))+
              geom_text(aes(x = low10_line,y=0.4,label = "Low 10% PRS"))+
              geom_density(alpha=0.4, position = 'identity')

  ggsave(dist_plot,file=tmp_dist_plot,height = 8,width  = 8)

  #dev.off()

########## Plotting Coefficients on Odds Scale ##########

  # Fit regression model
  prs_glm <- glm(Pheno ~ percentile,data = PRS.pheno.plot,family = 'binomial')

  # Put results in data.frame
  summs <- prs_glm %>% summary()

  # Get point estimates and SEs
  results <- bind_cols(coef(prs_glm),summs$coefficients[, 2]) %>%
  setNames(c("estimate", "se"))  %>%
  mutate(percentile = 1:10)

  # Your coefficients are on the log odds scale, such that a coefficient is
  # log(odds_y==1 / odds_y == 0). We can exponentiate to get odds instead.
  results_odds <- results %>% mutate(across(.cols = -percentile, ~ exp(.x)))

  # Need SEs on the odds scale too
  results_odds <- results_odds %>%
    mutate(var_diag = diag(vcov(prs_glm)),
           se = sqrt(estimate ^ 2 * var_diag))

  out_df_or<-results_odds
  out_df_or<-as.data.frame(out_df_or)

  tmp_out_or<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],".OR.data.txt")
  fwrite(out_df_or,tmp_out_or,sep="\t",col.names = T)


  # Plot with +/- 1 SE

  tmp_or_plot<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.OR-Plot.png',collapse = '')

  or_plot =ggplot(results_odds, aes(x = as.factor(percentile), y = estimate, color)) +
                    geom_point(stat = "identity", size=3,color = "black") +
                    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
                    geom_errorbar(aes(ymin = estimate - se, ymax = estimate + se), width = 0.4) +
                    ggtitle("Odds Ratio in 1st to 10th PRS score") +
                    xlab("PRS Quantile") +
                    ylab("Odds")
  ggsave(or_plot,file=tmp_or_plot,height = 8,width  = 8)

########## Quantiles plot ##########

  tmp_quant_plot<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.quantiles.png',collapse = '')

 quant_plot<- PRS.pheno.plot %>%
   count(percentile, Pheno) %>%       
   group_by(percentile) %>%
   mutate(pct= prop.table(n) * 100) %>%
   ggplot() + aes(percentile, pct, fill=Pheno) +
                  geom_bar(stat="identity") +
                  xlab("Quantiles") +
                  ylab("Ratio of Case/Control") +
                  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
                  position=position_stack(vjust=0.5)) +
                  ggtitle("PRS in 1st to 10th distribution") +
                  theme_bw()
    ggsave(quant_plot,file=tmp_quant_plot,height = 8,width  = 8)

########## Two sample T test: Case and control ########## 

  case<-subset(PRS.pheno.plot,PRS.pheno.plot$Pheno=="Case")
  ctrl<-subset(PRS.pheno.plot,PRS.pheno.plot$Pheno=="Control")

# Two sample T-test
  t.test.result<-tidy(t.test(case$SCORE,ctrl$SCORE))
  cat('\nP-value of PRS Case/Control distribution T-test is:',t.test.result$p.value)
  cat('\n')
  cat('\n')
  tmp_t_test<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.T-test.txt',collapse = '')
  fwrite(t.test.result,tmp_t_test,sep="\t",col.names = T)

# Wilcoxon Whitney U test , specify alternative="less"
  u.test.result<-tidy(wilcox.test(case$SCORE,ctrl$SCORE,alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE))
  cat('\nP-value of PRS Case/Control distribution U-test is:',u.test.result$p.value)
  cat('\n')
  cat('\n')
  tmp_u_test<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.U-test.txt',collapse = '')
  fwrite(u.test.result,tmp_u_test,sep="\t",col.names = T)

############################### Ready GLM PRS Model Data (Drop NA's) ########################################

  PRS.pheno.model=PRS.result.model
  PRS.pheno.model=na.omit(PRS.pheno.model)

  colnames(PRS.pheno.model)<-c("PatientID","IID","Sex","Age","Pheno","SCORE","PC1","PC2","PC3","PC4")
  
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

  tmp.tarin.data<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.AUC.train-data.txt',collapse = '')
  fwrite(PRS.train.pheno,tmp.tarin.data,sep="\t",col.names=T)

  tmp.test.data<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.AUC.test-data.txt',collapse = '')
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

  mod_1 <- glm( Pheno ~ SCORE, data=PRS.train.pheno, family="binomial")
  mod_2 <- glm( Pheno ~ SCORE+Age+Sex, data=PRS.train.pheno, family="binomial")
  mod_3 <- glm( Pheno ~ SCORE+Age+Sex+PC1+PC2+PC3+PC4, data=PRS.train.pheno, family="binomial")

##################################### Model 1 (Base only) #####################################

  library(pROC)

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

  out_final<-rbind(out_final_1,out_final_2,out_final_3)
  
  out_final$PGSID <- PGS_NUM[j]

  out_final<-out_final[,c(9,8,7,1:6)]

  tmp_out<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.Performance.txt')
  
  fwrite(out_final,tmp_out,sep="\t",col.names = T)

  ### Model plot ###

    tmp_model_plot<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.AUCs.png',collapse = '')

    png(tmp_model_plot,height = 500,width  = 500)

    test_1_roc = roc(PRS.test.pheno$Pheno ~ test_1_prob)
    test_2_roc = roc(PRS.test.pheno$Pheno ~ test_2_prob)
    test_3_roc = roc(PRS.test.pheno$Pheno ~ test_3_prob)

    plot(test_1_roc,print.auc = TRUE, print.auc.y = .4)
    plot(test_2_roc,print.auc = TRUE, print.auc.y = .3 , col='red' ,add=TRUE)
    plot(test_3_roc,print.auc = TRUE, print.auc.y = .2 , col='blue',add=TRUE)
    
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

################# Table 1 of PRS Sample in GLM model (Dropped NAs) #################

  cat('\n')
  cat('\n##### Table 1 of PRS Sample in GLM model (Dropped NAs) #####')
  cat('\n')

  table_PRS=PRS.pheno.model
  #table_PRS$Pheno<-recode_factor(table_PRS$Pheno,"1"="Control","2"="Case")
  #table_PRS$Sex<-recode_factor(table_PRS$Sex,"1"="Male","2"="Female")

  cov = c('Pheno','Sex','Age')
  table_PRS_print <- table_PRS %>% select(cov)

  tmp_PRS_table1<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.table1.html',collapse = '')

  table1 <- 
    tbl_summary(
      table_PRS_print,
      by = Pheno, # split table by group
      missing = "no" # don't list missing data separately
    ) %>%
    add_n() %>% # add column with total number of non-missing observations
    add_p() %>% # test for a difference between groups
    modify_header(label = "**Covarites**") %>% # update the column header
    bold_labels() 

  table1%>%
    as_gt() %>%
    gt::gtsave(filename = tmp_PRS_table1)
    
  cat('\n')
  cat('\nPRS Pheno table1 is in:',tmp_PRS_table1)
  cat('\n')
  cat('\n')

######################## Prevalence plot and Table ################################

  Prevalence.plot<-PRS.pheno.plot
  Prevalence.plot$Pheno<-recode_factor(Prevalence.plot$Pheno,"Control"="0","Case"="1")
  Prevalence.plot$Pheno<-as.numeric(as.character(Prevalence.plot$Pheno))
  Prevalence.plot$percentile<-as.factor(Prevalence.plot$percentile)
  Prevalence.plot$PrevalenceGroup<-ntile(Prevalence.plot$SCORE,100)

  PrevalencePlot_data<-Prevalence.plot %>% 
  group_by(PrevalenceGroup) %>% 
  summarise(Prevalence = sum(Pheno)/n())

  tmp_PrevalencePlot_data<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.Prevalence-data.txt')
  fwrite(PrevalencePlot_data,tmp_PrevalencePlot_data,sep="\t",col.names = T)


  tmp_Prevalence_plot<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.Prevalence.png')
        
  Prevalence_plot<-ggplot(PrevalencePlot_data, aes(x=PrevalenceGroup, y=Prevalence)) + 
  labs( x = "Percentile of PRS", y = "Prevalence",title ="Prevalence of Disease")+
  geom_point()

  ggsave(Prevalence_plot,file=tmp_Prevalence_plot,height = 8,width  = 8)


######################## Cumulative Risk plot and Table ################################

    Cumulative.plot<-PRS.pheno.plot


    Cumulative.plot$Pheno<-recode_factor(Cumulative.plot$Pheno,"Control"="0","Case"="1")
    Cumulative.plot$Pheno<-as.numeric(as.character(Cumulative.plot$Pheno))
    Cumulative.plot$percentile<-as.factor(Cumulative.plot$percentile)
    Cumulative.plot$PrevalenceGroup<-ntile(Cumulative.plot$SCORE,100)

    Cumulative.plot<-Cumulative.plot[order(Cumulative.plot$SCORE),]
    Cumulative.plot$Index<-1:nrow(Cumulative.plot)
    cumplot_glm <- fitSmoothHazard(Pheno ~ Age+SCORE,data = Cumulative.plot, time = "Age", ratio = 10)

    #summary(cumplot_glm)

    group_25<-subset(Cumulative.plot,Cumulative.plot$PrevalenceGroup=="25")
    group_50<-subset(Cumulative.plot,Cumulative.plot$PrevalenceGroup=="50")
    group_75<-subset(Cumulative.plot,Cumulative.plot$PrevalenceGroup=="75")
    group_100<-subset(Cumulative.plot,Cumulative.plot$PrevalenceGroup=="100")

    sample_25<-min(group_25$Index)
    sample_50<-min(group_50$Index)
    sample_75<-min(group_75$Index)
    sample_100<-min(group_100$Index)

    smooth_risk_model <- absoluteRisk(object = cumplot_glm, newdata = Cumulative.plot[c(sample_25,sample_50,sample_75,sample_100),])

    risk_data<-as.data.frame(smooth_risk_model)
    colnames(risk_data)<-c("Age-Time","PRS-25%","PRS-50%","PRS-75%","PRS-100%")

    tmp_risk_data<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.Cumulative_Risk-data.txt')
    fwrite(risk_data,tmp_risk_data,sep="\t",col.names = T,row.names=F)

    tmp_risk_plot<-paste0('./output/',Pheno_col,'/',PGS_NUM[j],'/plot/',PGS_NUM[j],'.Cumulative_Risk.png')

    smooth_risk <- absoluteRisk(object = cumplot_glm, newdata = Cumulative.plot[c(sample_25,sample_50,sample_75,sample_100),])

    png(tmp_risk_plot,height = 800,width  = 800)
    fullplot<-plot(smooth_risk,
          id.names = c("PRS 25%","PRS 50%","PRS 75%","PRS 100%"), 
          legend.title = "PRS", 
          xlab = "Age (Years Old)", 
          ylab = "Cumulative Incidence (%)")
    print(fullplot)
    dev.off()

 }
    
}

