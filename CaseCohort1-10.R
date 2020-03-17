CaseCohort <-function(dat, dat.id, event, id , y){
  #set.seed(24)
  datcc.id <- dat.id
  datcc.id$sample <- rbinom(nrow(datcc.id), 1, (0.1))
  datcc.id$sample2 <- ifelse(datcc.id[[event]] == 1, 1,
                             ifelse(datcc.id$sample == 1 , 1, 0))

  # Merge this also into the dataset in the long format
  test.id <-  datcc.id[c(id, "sample2")]
  test <- merge(test.id, dat)
  datcc <- test[order(test$id),]
  
  # Data for Fit II
  datcc$y2 <- ifelse(datcc$sample2 == 0 , NA, datcc[[y]] )
  # Data for Fit III
  datcc3 <- subset(datcc, datcc$sample2 == 1)
  datcc3.id <- subset(datcc.id, datcc.id$sample2 == 1)

  return(list(datcc, datcc3, datcc3.id))

}

