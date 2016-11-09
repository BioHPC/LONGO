LONGOquotient <- function(
                           DATA.DF,
                           JS_DATA.DF,
                           CONTROL_COLUMN_INDEX
                         ){

  max_js <- max(JS_DATA.DF[nrow(JS_DATA.DF),2:ncol(JS_DATA.DF)])
  for (i in 2:ncol(JS_DATA.DF)) {
    JS_DATA.DF[,i] <- JS_DATA.DF[,i]/max_js
  }

  JS_DATA.DF[,CONTROL_COLUMN_INDEX] <- 0
  LQ_data <- JS_DATA.DF[!(JS_DATA.DF[,1]<150),]
  LQ_sign_data <- DATA.DF[!(JS_DATA.DF[,1]<150),]
  LQ_sign <- data.frame(matrix(0, nrow = 1, ncol = ncol(DATA.DF)))

  colnames(LQ_sign) <- colnames(LQ_data)
  for (i in 2:ncol(LQ_sign_data)) {
    avg <- mean(LQ_sign_data[,i])

    avg_control <- mean(LQ_sign_data[,CONTROL_COLUMN_INDEX])
    if (avg >= avg_control) {
      LQ_sign[1,i] <- 1
    } else {
    LQ_sign[1,i] <- -1
    }
  }
  LQ_sign <- LQ_sign[,-1]
  for (i in 2:ncol(LQ_data)) {
    min_LQ <- min(LQ_data[,i])
    for (j in 1:nrow(LQ_data)) {
      LQ_data[j,i] <- LQ_data[j,i] - min_LQ
    }
  }
  z <- data.frame(matrix(0, nrow = 1, ncol = ncol(LQ_data)))
  for (i in 2:ncol(LQ_data)) {
    z[1,i] <- LQ_data[nrow(LQ_data),i]-LQ_data[1,i]
  }
#  z <- as.matrix(z)
 # z<<-z
  colnames(z) <- colnames(DATA.DF)[2:ncol(DATA.DF)]
#  z <- as.data.frame(z)
  z <- z[,-1]
#  labels <- colnames(z)
#  labels <- as.matrix(labels)
  for (i in 1:ncol(z)) {
    z[1,i] <- z[1,i]*LQ_sign[1,i]
  }
  z <- t(z)
  z <- z[order(z[,1]),]
  z <- t(z)
  return(z)
}



