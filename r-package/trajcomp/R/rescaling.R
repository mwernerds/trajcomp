# input: trajectory (x,y[,t]),
#         rescaling length
# output: rescaled time series


rescale <- function(data, n){
  
  library(TSclust)

  nOld = nrow(data);

  # check whether there is a time column
  if("t" %in% colnames(data)){
    print("time column exists")
  }
  # else add time column (required for linear interpolation): 1,2,..., nOld
  else{
    data = cbind(data,seq(1,nOld))
  }
 
  dataRescaled <- data.frame();

  # time series has the required length
  if(nOld == n){
    dataRescaled <- data;
  }

  # time series is too short -> piecewise linear interpolation (PLI)
  else if(nOld < n){
    dataRescaled <- doPLI(data, n);

  }

  # time series is too long -> Piecwise Aggregate Approximation (PAA)
  else if(nOld > n){
    dataRescaled <- TSclust::PAA(data[,1], n);
    dataRescaled <- cbind(dataRescaled, TSclust::PAA(data[,2], n));
    dataRescaled <- cbind(dataRescaled, TSclust::PAA(data[,3], n));
  }
  
  return(dataRescaled);
}
