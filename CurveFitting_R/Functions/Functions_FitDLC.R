## functions to fit DLC

funfit <- function(params, t1) {
  t2 <- (1 / (1 + exp((params[1] - t1) / params[2]))) - (1 / (1 + exp((params[3] - t1) / params[4])))
  return(t2)
}

objective_function <- function(params, ndvi, t1) {
  predicted <- funfit(params, t1)
  sum((ndvi - predicted) ^ 2)
}

funfitspring <- function(x){
  if(values=="daily"){
    time=1:365
    t2 <- 1/(1+exp((x[1]-time)/x[2]))
  }else{
    time=seq(1,365,by=7)
    t2 <- 1/(1+exp((x[1]-time)/x[2]))
  }
  return(t2)
}
funIRG <- function(x){
  if(values=="daily"){
    t <- 1:365
  }else{
    t <- seq(1,365,by=7)
  }
  t2 <- c((diff(x)/diff(t)),NA)
  t2[t2 < 0] <- 0
  if(sclIRG==TRUE){
    if(all(is.na(t2))==FALSE){
      t2 <- (t2-min(t2, na.rm=TRUE))/(max(t2, na.rm=TRUE)-min(t2, na.rm=TRUE))
    }
  }
  return(t2)
}


dlc_coefs <-function(ndvi_data){
  
  model <- y ~ (1/(1+exp((xmidS-time)/scalS))) - (1/(1+exp((xmidA-time)/scalA)))
  lowers <- c(Xmid1=30, Scale1=1, Xmid2=167, Scale2=1)
  uppers <- c(Xmid1=228, Scale1=30, Xmid2=288, Scale2=35)
  time <- as.numeric(strftime(ndvi_data$date, format="%j", tz="GMT"))
  
  res <- try(nls(model, data=data.frame(y=ndvi_data$rescale, time=time),
                 start=list(xmidS=120, scalS=20, xmidA=240, scalA=20),
                 lower=c(xmidS=1, scalS=1, xmidA=200, scalA=1),
                 upper=c(xmidS=365, scalS=100, xmidA=365, scalA=100),
                 control=nls.control(maxiter=100, minFactor=1/1024),
                 algorithm="port"), silent=TRUE)

  if(class(res)=="try-error"){
    
    return.list <- list(predicted_values = NULL, dlcparams = NULL, irg = NULL)
    
  }else{

    s <- try(summary(res), silent=T)
    dlcparams <- c(summary(res)$coefficients[,1],summary(res)$coefficients[,2])
    
    dailytime <- 1:365
    
    ## predicted for the year
    temp <- data.frame(predict_Vals = (1/(1+exp((dlcparams[1]-dailytime)/dlcparams[2])))-(1/(1+exp((dlcparams[3]-dailytime)/dlcparams[4]))))
    temp$jul <- 1:365

    ## predicted for the spring
    temp_spring <- 1/(1+exp((dlcparams[1]-dailytime)/dlcparams[2]))
    
    ## IRG
    irg <- c((diff(temp_spring)/diff(dailytime)),NA)
    irg[irg < 0] <- 0
    
    
    ## get the unscaled predictions from the DLC
    min_ndvi <- min(ndvi_data$NDVI)
    max_ndvi <- max(ndvi_data$NDVI)
    predicted_values <- min_ndvi + (temp$predict_Vals - min(temp$predict_Vals)) * (max_ndvi - min_ndvi) / (max(temp$predict_Vals) - min(temp$predict_Vals))

    return.list <- list(predicted_values = predicted_values, dlcparams = dlcparams, irg = irg)
    
  }
  
  
  return(return.list)
  
  
}

