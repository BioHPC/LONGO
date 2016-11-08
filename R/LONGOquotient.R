#LONGOquotient <- function(
#                            DATA.DF,
#                            JS_DATA.DF
#                            #CONTROL COLUMN NEEDED?
 #                         ){

#data.df = simp.df
#y=jsdistance

#max_y <- max(y[nrow(y),2:ncol(y)])
#for (i in 2:ncol(y)) {
#  y[,i] <- y[,i]/max_y
#  #y[,i] <- y[,i]/max(y[,i])
#}

#y[,control_index] <- 0

#LQ_data <- y[!(y[,1]<150),]
#LQ_sign_data <- data.df[!(y[,1]<150),]
#LQ_sign <- data.frame(matrix(0, nrow = 1, ncol = ncol(data.df)))

#colnames(LQ_sign) <- colnames(LQ_data)
#for (i in 1:ncol(LQ_sign_data)) {
#  avg <- mean(LQ_sign_data[,i])
#  avg_control <- mean(LQ_sign_data[,control_index])
#  if (avg >= avg_control) {
#    LQ_sign[1,i] <- 1
#  } else {
#    LQ_sign[1,i] <- -1
#  }
#}
#LQ_sign <- LQ_sign[,-1]
#for (i in 2:ncol(LQ_data)) {
#  min_LQ <- min(LQ_data[,i])
#  for (j in 1:nrow(LQ_data)) {
#    LQ_data[j,i] <- LQ_data[j,i] - min_LQ
#  }
#}



#plot. 1



#cutoff


#}
