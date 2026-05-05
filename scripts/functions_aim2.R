#Classify transitions between nutrititonal statuses
classify_transitions <- function(stunted_list) {
  unique_stunted <- unique(stunted_list)
  if (length(unique_stunted) == 1) {
    return(unique_stunted)
  } else {return(paste(unique_stunted,collapse=" , ")) }
}


#Calculates variance and bias of a RF model
bias_variance = function(dataset,iterations,important_species_df,non_important_species_df){
  results_bias_variance = data.frame()
  index_train = sample(nrow(dataset),nrow(dataset)*0.7,replace = F)
  dataset_test = dataset[-index_train,]
  dataset_train = dataset[index_train,]
  M = iterations
  step = ncol(dataset_train) %/% 40
  length_important = length(important_species_df$species)
  length_non_important = length(non_important_species_df$species)
  
  #Features (Variables) in the dataset 
  for(i in seq(from=1,to=length(c(important_species_df$species,non_important_species_df$species)),by=step)){
    set.seed(50 + i)
    
    #predictions based on the current number of features (i)
    test_predictions_matrix = matrix(NA, nrow=nrow(dataset_test), ncol=M)
    
    #iterations of RF. Each iteration is a new sample of the dataset 
    for(m in 1:M){
      
      #Get trainining dataset with desired number of features
      if(i > length_important){
        sample_non_important_totrain=sample(non_important_species_df$species,i-length_important,replace = F)
        train = dataset_train[,names(dataset_train) %in% c(important_species_df$species[1:i],sample_non_important_totrain,"age")]
      }else{
        train = dataset_train[,names(dataset_train) %in% c(important_species_df$species[1:i],"age")]
      }
      
      
      #Random Forest Implementation 
      rf_training=randomForest(y=train$age,x=train[,-c(ncol(train))],subset=1:nrow(train))
      
      #Get test dataset 
      feature_names = names(train)[-ncol(train)]
      test_healthy = predict(rf_training, dataset_test[, feature_names])
      #Update predictions matrix for i number of features
      test_predictions_matrix[, m] <- test_healthy
    }
    
    #Bias: Bias^2=(E[f^hat(x)]−f(x))^2 ; E[f_hat(x)]: average predictions over multiple training sets
    #Variance: E[(f_hat(x)-E[f_hat(x)])^2]
    
    #Real data f(x)
    true_age = dataset_test$age
    
    # Average Prediction over multiple training sets (E[f_hat(x)])
    avg_pred = rowMeans(test_predictions_matrix)
    
    #Bias squared
    bias_squared = mean((avg_pred - true_age)^2)
    
    #Variance is calculated using the mean of var() for the variance of the predictions from all the models (rows)
    variance = mean(apply(test_predictions_matrix, 1, var))
    
    #Store bias and variance for all number of features
    results_bias_variance = rbind(results_bias_variance, data.frame(num_features=length(feature_names),bias_sq=bias_squared,variance = variance))
  }
  
  return(results_bias_variance)
}

#Runs CV  
cross_val_predictors = function(dataset,important_species_df,non_important_species_df,number_background_pred,number_iterations,cv_fold){
  sample_non_important=sample(non_important_species_df$species,number_background_pred,replace = F)
  result <- rfcv(dataset[,names(dataset) %in% c(important_species_df$species,sample_non_important)], dataset$age, cv.fold=cv_fold)
  n_steps <- length(result$n.var)
  error_matrix <- matrix(NA, nrow=n_repeats, ncol=n_steps)
  for (i in 1:n_repeats) {
    set.seed(50 + i)
    result <- rfcv(healthy_matrix_train[,names(healthy_matrix_train) %in% c(important_species$species,sample_non_important)], healthy_matrix_train$age, cv.fold=10) 
    error_matrix[i, ] <- result$error.cv
  }
  
  error_mean <- colMeans(error_matrix)
  error_sd <- apply(error_matrix, 2, sd)
  n_vars <- result$n.var
  df <- data.frame(n_vars=n_vars, error_mean=error_mean, error_sd=error_sd)
  
  return(list(error_matrix,df))
}


#Select train and test datasets for a matrix. 
selectDatasets=function(matrix,infants_to_train,predictors){
  set.seed(1208)
  
  trainIndex <- which(matrix$subjectID.x %in% infants_to_train)
  trainData <- matrix[trainIndex, ]
  testData <- matrix[-trainIndex, ]
  all_testData <- matrix[!(matrix$subjectID.x %in% infants_to_train),]
    
  trainData<-trainData[,colnames(trainData) %in% predictors]
  testData <-testData[,colnames(testData) %in% predictors]

  trainData$sample=matrix[trainIndex,"sample"]
  testData$sample=matrix[-trainIndex,"sample"]
  
  testData = testData %>% unnest(sample) %>% mutate(sample = as.character(sample))
  trainData = trainData %>% unnest(sample) %>% mutate(sample = as.character(sample))
  
  trainData = trainData %>% left_join(all_children[,c(3,7)],
                                      by = c("sample" = "accession"))
  testData = testData %>% left_join(all_children[,c(3,7)],
                                    by = c("sample" = "accession"))
  return(list("train_healthy"=trainData,
              "test_healthy"=testData,
              "test_all"=all_testData))
}


#Select datasets for vOTUs.(taking into account shannon and replication lifestyle)
#To select train and test
selectDatasets_1=function(vOTU_matrix,infants_to_train,lifestyle,shannon=F,lifestyle_df){
  set.seed(1208)
  healthy_df = vOTU_matrix[vOTU_matrix$trajectory=="Consistent normal growth",]
   
  trainIndex <- which(healthy_df$subjectID.x %in% infants_to_train)
  trainData <- healthy_df[trainIndex, ]
  testData <- healthy_df[-trainIndex, ]
  
  #Define which contigs go into the dataset based on lifestyle (vector)
  list_contigs = unique(lifestyle_df$contig[lifestyle_df$replication_cycle %in% lifestyle])
  trainData<-trainData[,colnames(trainData) %in% list_contigs]
  testData <-testData[,colnames(testData) %in% list_contigs]
  
  trainData$sample=healthy_df[trainIndex,"sample"]
  testData$sample=healthy_df[-trainIndex,"sample"]
  
  testData = testData %>% unnest(sample) %>% mutate(sample = as.character(sample))
  trainData = trainData %>% unnest(sample) %>% mutate(sample = as.character(sample))
  
  trainData = trainData %>% left_join(all_children[,c(3,7)],
                                      by = c("sample" = "accession"))
  testData = testData %>% left_join(all_children[,c(3,7)],
                                    by = c("sample" = "accession"))
  
  #Define if shannon makes part of the dataset
  if(shannon){
    if(length(lifestyle)==4){
      data=df_alpha_all
      trainData = trainData %>% left_join(data[,c(1,3)],by="sample")
      testData = testData %>% left_join(data[,c(1,3)],by="sample")
    }else{
      data = df_alpha_all_lifeCycle[df_alpha_all_lifeCycle$replication_cycle %in% lifestyle,]
      trainData = trainData %>% left_join(data[,c(1,4)],by="sample")
      testData = testData %>% left_join(data[,c(1,4)],by="sample")
    }
  }
  return(list(trainData,testData))
}

#Select datasets for other datasets (ex. microdiversity) IAM NOT USING THIS 
selectDatasets_other=function(dataset){
  set.seed(1208)
  healthy_df = dataset[dataset$stunted_status=="Healthy",]
  infants_healthy = healthy_df$subjectID.x
  sampled_infants <- sample(unique(infants_healthy),
                            size = length(unique(infants_healthy))*0.8)
  
  trainIndex <- which(healthy_df$subjectID.x %in% sampled_infants)
  trainData <- healthy_df[trainIndex, ]
  testData <- healthy_df[-trainIndex, ]
  
  return(list(trainData,testData))
}

# Create a reference grid for EMMs across age ONLY FOR vOTUs
emmeans_alpha = function(gam_model,age_intervals,data_gam){
  ref_grid_richness <- ref_grid(gam_model, at = list(age = age_intervals, total_reads = mean(data_gam$total_reads)))
  emm <- emmeans(ref_grid_richness, ~ trajectory | age)
  emm_df <- as.data.frame(emm)
  return(emm_df)
}
emmeans_shannon_across = function(gam_model,age_intervals,data_gam){
  ref_grid_richness <- ref_grid(gam_model, at = list(age = age_intervals, total_reads = mean(data_gam$total_reads)))
  emm <- emmeans(ref_grid_richness, ~ trajectory | age)
  emm_df <- as.data.frame(emm)
  return(list(emm_df,emm))
}

analyze_age_pairwise <- function(model, data, age_range, type_label = "data") {
  emm_out <- emmeans_shannon_across(model, age_range, data)
  pairwise_obj <- pairs(emm_out[[2]], by = "age", adjust = "tukey")
  
  
  pairwise_df <- summary(pairwise_obj, infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    tidyr::separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
    dplyr::relocate(age, group1, group2) %>%
    dplyr::mutate(
      significant = p.value < 0.05,
      source = type_label  # Adds a column to identify if it's 'p' or 'b'
    )
  summary_stats <- pairwise_df %>%
    dplyr::group_by(age) %>%
    dplyr::summarise(
      n_contrasts = n(),
      min_p_adj   = min(p.value, na.rm = TRUE),
      max_abs_diff = max(abs(estimate), na.rm = TRUE),
      any_sig = any(p.value < 0.05, na.rm = TRUE),
      source = type_label
    ) %>%  dplyr::arrange(age)
  
  return(details = pairwise_df)
}

#Run Boruta and get important and non important taxa for maturation
feature_selection = function(train_df,predictors,tree_number=500,iterations=500){
  #Run Boruta
  Boruta.age <- Boruta(train_df$age~.,data=train_df[,names(train_df) %in% predictors],doTrace=2,ntree=tree_number,maxRuns = iterations)
  #GET IMPORTANT SPECIES
  important_species = as.data.frame(Boruta.age$ImpHistory) %>% pivot_longer(cols = everything(),names_to = "species",values_to = "importance")  %>%  group_by(species) %>% summarise(mean_importance=mean(importance)) %>% left_join( as.data.frame(Boruta.age$finalDecision) %>% mutate(species = rownames(as.data.frame(Boruta.age$finalDecision))), by = "species") %>% filter(`Boruta.age$finalDecision` != "Rejected") %>% select(species)
  #GET IMPORTANT SPECIES
  non_important_species = as.data.frame(Boruta.age$ImpHistory) %>% pivot_longer(cols = everything(),names_to = "species",values_to = "importance")  %>%  group_by(species) %>%summarise(mean_importance=mean(importance)) %>% left_join( as.data.frame(Boruta.age$finalDecision) %>% mutate(species = rownames(as.data.frame(Boruta.age$finalDecision))), by = "species") %>% filter(`Boruta.age$finalDecision` == "Rejected") %>% select(species)
  
  
  #Importance and correlation with age 
  R2_df= train_df %>% pivot_longer(cols=predictors,names_to="species",values_to = "rel.abundance") %>% select(species,rel.abundance,age) %>% group_by(species) %>% summarize(cor_withage= sign(cor(rel.abundance,age)),r2=cor(rel.abundance,age)^2,sign_R2=r2*cor_withage)
  
  df_importance = as.data.frame(Boruta.age$ImpHistory) %>% pivot_longer(cols = everything(),names_to = "species",values_to = "importance")  %>%  group_by(species) %>% summarize(mean_importance=mean(importance), sd_importance = sd(importance)) %>% left_join(as.data.frame(Boruta.age$finalDecision) %>% mutate(species = rownames(as.data.frame(Boruta.age$finalDecision))), by = "species") %>% left_join(R2_df,by="species") %>% filter(`Boruta.age$finalDecision` != "Rejected") 
  
  return(df_importance)
}

#PLot importance and R2 with age for predictors selected by Boruta.
plots_maturation = function(df_importance,train_df){
  
  #R2 and importance index plot
  max_abs_R2 <- max(abs(df_importance$sign_R2), na.rm = TRUE)
  min_abs_R2 <- min(abs(df_importance$sign_R2), na.rm = TRUE)
  
  symm_limits <- c(-0.5, max_abs_R2)
  neg_colors <- viridis::mako(100,begin = 0.3)
  pos_colors <- viridis::magma(100,begin = 0.6)
  all_colors <- c(neg_colors,rev(pos_colors))
  #barplot(rep(1,length(all_colors)), col=all_colors)
  
  
  pred_imp_plot = df_importance %>% ggplot(aes(x=reorder(species,mean_importance),y=mean_importance,color=sign_R2))+ geom_point() +  geom_line() + geom_errorbar(aes(ymin=mean_importance - sd_importance, ymax=mean_importance + sd_importance), width=0.2) + coord_flip() + ylab("Mean Decrease Gini impurity") + xlab("")+theme_classic2() + scale_color_gradientn(colors = all_colors,limits = symm_limits,values = scales::rescale(c(symm_limits[1], 0, symm_limits[2])),name=expression("sign *  R"^2))+theme(strip.placement  = "outside",strip.text.y = element_text(angle=0))
  
  #Abundance maturation patterns
  abundance_moreimportant = train_df[,names(train_df) %in% df_importance$species]
  abundance_moreimportant = cbind(abundance_moreimportant,train_df[,c("sample","age")])
  abundance_moreimportant_melt= reshape2::melt(abundance_moreimportant,id.vars = c("sample","age"))
  
  abundance_moreimportant_melt$sample = factor(abundance_moreimportant_melt$sample)
  abundance_moreimportant_melt$variable = factor(abundance_moreimportant_melt$variable)
  abundance_moreimportant_melt$value = as.numeric(as.character(abundance_moreimportant_melt$value))
  
  abundance_moreimportant_melt = abundance_moreimportant_melt %>% group_by(age,variable) %>% summarise(mean_relab=median(value))
  abundance_moreimportant_melt$age=round(abundance_moreimportant_melt$age)
  
  abundance_important_age_matrix <- reshape2::dcast(abundance_moreimportant_melt, variable~age) 
  
  rownames(abundance_important_age_matrix)<-abundance_important_age_matrix$variable
  abundance_important_age_matrix<-data.matrix(abundance_important_age_matrix[,-1])
  abundance_important_age_matrix[is.na(abundance_important_age_matrix)] <- 0 
  dim(abundance_important_age_matrix)
  
  abundance_important_age_matrix <- abundance_important_age_matrix[
    apply(abundance_important_age_matrix, 1, function(x) sd(x) > 0), ]
  
  ages <- as.numeric(colnames(abundance_important_age_matrix))
  age_bins <- cut(ages, breaks = seq(0, max(ages), by = 50), right = FALSE)
  
  collapsed <- sapply(levels(age_bins), function(bin) {
    cols_in_bin <- which(age_bins == bin)
    rowSums(abundance_important_age_matrix[, cols_in_bin, drop=FALSE], na.rm=TRUE)
  })
  
  rownames(collapsed) <- rownames(abundance_important_age_matrix)
  
  collapsed[is.na(collapsed)] <- 0 
  
  collapsed <- collapsed[apply(collapsed, 1, function(x) sd(x) > 0), ]
  
  plot_abundance_matrix =pheatmap(as.matrix(collapsed),color=viridis::magma(400,begin=-0,end=1,direction = 1),scale = "row",cluster_cols = F,fontsize = 20,clustering_method="average",clustering_distance_rows="euclidean")

  return(list(pred_imp_plot,plot_abundance_matrix,collapsed))
}

normal_growth_maturation = function(train_df,predictors,test_df_allTrajectories,organism){
  #### Healthy maturation model #### 
  tic()
  rf_training=randomForest(y=train_df$age,x=train_df[,names(train_df) %in% predictors],subset=1:nrow(train_df))
  toc()
  
  healthy_df_test = test_df_allTrajectories[test_df_allTrajectories$trajectory=="Consistent normal growth",]
  #hist(as.numeric(healthy_matrix_test$age))
  
  test_healthy=predict(rf_training,healthy_df_test[,names(healthy_df_test) %in% predictors])
  
  predictions<-data.frame("sample"=healthy_df_test[,c("sample")],"pred_age"=test_healthy) %>% left_join(all_children,by=c("sample"="accession"))
  
  #ALL STATUSES - RANDOM FOREST
  
  test_all=predict(rf_training,test_df_allTrajectories[,names(test_df_allTrajectories) %in% predictors])
  
  predictions<-data.frame("sample"=test_df_allTrajectories[,c("sample")],"pred_age"=test_all) %>% left_join(all_children,by=c("sample"="accession"))
  
  square_error_df=NULL
  for(i in unique(predictions$trajectory)){
    status = i 
    prediction_data = predictions[which(predictions$trajectory == i),]
    mae <- mean(abs(prediction_data$pred_age - prediction_data$age))
    mse <- mean((prediction_data$pred_age - prediction_data$age)^2)
    rmse <- sqrt(mse)
    square_error = data.frame("trajectory"=i,"squared_error"=sqrt((prediction_data$pred_age - prediction_data$age)^2),"organism"=organism,
                              rmse,mse,mae)
    square_error_df=rbind(square_error_df,square_error)
    result = paste("RMSE_",i," = ",rmse,sep="")
    print(result)
  }
  
  
  library(mgcv)
  predictions_healthy = predictions %>% filter(trajectory=="Consistent normal growth")
  predictions_healthy$subjectID.x=factor(predictions_healthy$subjectID.x)
  #gam_model <- gam(pred_age ~  s(age,k=5) + s(subjectID.x, bs = "re"), data = predictions_healthy, family = Gamma(link = "log"),gamma=1.5)
  #glmm_model <- lme4::lmer(pred_age ~ age + (1 | subjectID.x), data = predictions_healthy,REML = T)
  glmm_model <- lm(pred_age ~ age, data = predictions_healthy,REML = T)
  
  print(paste("AIC_glmm=",AIC(glmm_model)))
  #print(paste("AIC_gam=",AIC(gam_model)))
  #if(AIC(glmm_model)<=AIC(gam_model)){gam_model=glmm_model}
  
  gam_model = glmm_model
  library(gratia)
  # Create a new data frame over a range of age values
  newdata <- data.frame(age = seq(min(predictions_healthy$age),
                                  max(predictions_healthy$age),
                                  length.out = 100))
  
  # For random effect to work, you need a valid subjectID.x (pick a real one)
  # OR set it to a factor level seen in the training data
  newdata$subjectID.x <- factor(predictions_healthy$subjectID.x[10],
                                levels = levels(predictions_healthy$subjectID.x))
  
  # Predict fitted values
  newdata$pred_age <- predict(gam_model, newdata,type = "response")
  newdata$organism=organism
  predictions_healthy$organism=organism
  
  mm <- model.matrix(terms(gam_model), newdata)
  newdata$se <- sqrt(diag(mm %*% tcrossprod(vcov(gam_model), mm)))
  newdata$upper <- newdata$pred_age + (1.96 * newdata$se)
  newdata$lower <- newdata$pred_age - (1.96 * newdata$se)
  
  residuals <- predictions_healthy$age - predictions_healthy$pred_age
  ss_res <- sum(residuals^2)
  ss_tot <- sum((predictions_healthy$age - mean(predictions_healthy$age))^2)
  r2_value <- 1 - (ss_res / ss_tot)
  r2_label <- paste0("R^2 == ", round(r2_value, 3))
  
  y_axis=ifelse(organism=="Bacteria","Bacteriome","Virome")
  
plot_maturation_bacteria_healthy = ggplot()+geom_ribbon(data = newdata, aes(x = age, ymin = lower, ymax = upper),fill = "grey", alpha = 0.4)+geom_line(data= newdata, aes(x = age, y = pred_age),color = colorblind_palette[1], size = 1.2) +
    labs(y = paste(y_axis,"Age",sep=" "), x = "Infant Chronological Age") + theme_classic(base_size = 14)+scale_color_manual(values=colorblind_palette,name="")+ geom_point(data=predictions_healthy,aes(x=age,y=pred_age,colour = trajectory))+ scale_x_continuous(expand = expansion(mult = 0.1))+geom_abline(slope = 1,intercept = 0)+annotate("text",x=12,y=500,label = r2_label,parse=T,size=5,hjust=0)+theme(legend.position = "none")
  
  predictions_healthy$weaned =ifelse(predictions_healthy$age<200,"F","T")
  predictions_healthy$weaned_pred =ifelse(predictions_healthy$pred_age<200,"F","T")
  predictions_healthy$consistent=ifelse(predictions_healthy$weaned == predictions_healthy$weaned_pred,"Correctly assigned","Incorrectly assigned")
  
  print(table(predictions_healthy$consistent)[1]/sum(table(predictions_healthy$consistent)))
  
  plot_maturation_bacteria_healthy_weaning = ggplot()+geom_line(data= newdata, aes(x = age, y = pred_age),color = colorblind_palette[1], size = 1.2) + labs(y = paste(y_axis,"Age",sep=" "), x = "Infant Chronological Age") +  theme_classic(base_size = 14)+geom_point(data=predictions_healthy,aes(x=age,y=pred_age,color= consistent))+theme_pubr()+geom_abline(slope = 1,linetype=2)+geom_vline(xintercept = 200)+scale_color_manual(values=c("#C5A059","#0D1B2A"),name="Weaning Period\nPrediction")+annotate("text",x=12,y=500,label = r2_label,parse=T,size=5,hjust=0)+theme(legend.position = "none")
  
  data_comp_pred_data <- data.frame(age = predictions$age)
  
  data_comp_pred_data$subjectID.x <- factor(predictions_healthy$subjectID.x[10],
                                            levels = levels(predictions_healthy$subjectID.x))
  
  # Predict fitted values
  data_comp_pred_data$predicted_healthy <- predict(gam_model, data_comp_pred_data,type = "response")
  
  data_comp_pred_data$trajectory = predictions$trajectory
  data_comp_pred_data$predicted_age_rf= predictions$pred_age
  
  data_comp_pred_data$error= data_comp_pred_data$predicted_age_rf - data_comp_pred_data$predicted_healthy
  
  low_bond = quantile(data_comp_pred_data$error[data_comp_pred_data$trajectory=="Consistent normal growth"],probs = c(0.025))
  upp_bond = quantile(data_comp_pred_data$error[data_comp_pred_data$trajectory=="Consistent normal growth"],probs = c(0.975))
  median_healthy = quantile(data_comp_pred_data$error[data_comp_pred_data$trajectory=="Consistent normal growth"],probs = c(0.5))
  
  title_MG=paste(organism, sep="")
  maturity_gap = ggplot(data_comp_pred_data, aes(y = error, x = trajectory, color = trajectory)) +
    annotate("rect", ymin = low_bond, ymax = upp_bond, xmin = -Inf, xmax = Inf, fill = colorblind_palette[2], alpha = 0.2) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.4) + 
    scale_color_manual(values = colorblind_palette) + 
    xlab("") + ylab("Predicted age if healthy - Real age")+ geom_hline(yintercept = median_healthy, linetype = "dashed") + 
    stat_compare_means(method = "wilcox.test", ref.group = "Consistent normal growth",label = "p.signif",label.y = 150,na.rm = T,hide.ns = T) + 
    labs(title=title_MG)+ theme_classic(base_size = 20) + scale_y_continuous(limits = c(-350, 350)) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          title = element_text(size = 20))+
    guides(fill = guide_legend(nrow = 1, label.position = "top",byrow = TRUE))
  maturity_gap
  print(wilcox.test(data_comp_pred_data$error[data_comp_pred_data$trajectory=="Recovery to normal growth"], 
                    data_comp_pred_data$error[data_comp_pred_data$trajectory=="Consistent normal growth"]))
  
  return(list(gam_model_pred=newdata,healhy_pred=predictions_healthy,plot_maturation=plot_maturation_bacteria_healthy,plot_error=square_error_df,gam_model=gam_model,plot_maturation_weaning=plot_maturation_bacteria_healthy_weaning,plot_maturity_gap=maturity_gap))
  
}

r2_label = function(predictions_data,organism){
  residuals <- predictions_data$age - predictions_data$pred_age
  ss_res <- sum(residuals^2)
  ss_tot <- sum((predictions_data$age - mean(predictions_data$age))^2)
  r2_value <- 1 - (ss_res / ss_tot)
  r2_label <- paste0("R[",organism,"] ^2 ==", round(r2_value, 3))
  return(r2_label)
} 
