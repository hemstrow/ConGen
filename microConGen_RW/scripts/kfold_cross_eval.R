
dat <- read_genepop("Data/N50L100g2003.genepop")


test_samples <- sample(x = nsamps(dat), # how many samples do we have?
                       size = nsamps(dat)*.2) # how many do we want to hold back?


dat_num <- format_snps(dat, 
                       output = "sn")

dat_num <- dat_num[,-1] # remove SNP metadata--just need genotypes


kfold <- function(genos, phenos, k,
                  num.trees,
                  mtry){
  
  # shuff <- sample(length(phenos), length(phenos), FALSE)
  # 
  # genos <- genos[shuff,]
  # phenos <- phenos[shuff]
  
  parts <- split(1:length(phenos), sort(rep(1:k, length.out = length(phenos))))
  phenos <- as.factor(phenos)
  
  res <- vector("list", k)
  
  for(i in 1:length(parts)){
    cat("Working on part ", i, "\n")
    # remove samples and transpose (individuals in rows)
    GR_test <- t(genos[,parts[[i]]])
    GR_train <- t(genos[,-parts[[i]]])
    
    GR_train <- as.data.frame(GR_train)
    GR_test <- as.data.frame(GR_test)
    GR_train$pop <- phenos[-parts[[i]]]
    
    rf <- ranger(data = GR_train, num.trees = num.trees, 
                 mtry = mtry, # minus one for the pop column we added
                 dependent.variable.name = "pop",
                 keep.inbag = TRUE) # tell it that we are predicting population
    
    test_res <- forestError::quantForestError(rf, 
                                          X.train = GR_train[,-ncol(GR_train)], 
                                          X.test = GR_test,
                                          Y.train = GR_train$pop)
    
    
    # compare visually
    res[[i]] <- data.table(predicted = as.character(test_res$pred), observed = phenos[parts[[i]]], mcr = test_res$mcr)
  }
  
  res <- data.table::rbindlist(res)
  # res <- res[shuff,]
  
  cat("Total prediction success: ", mean(res$predicted == res$observed)*100)
  
  return(res)
}


kf <- kfold(dat_num, meta$genepop_pop, 10, num.trees = 100000, mtry = ncol(dat_num))
