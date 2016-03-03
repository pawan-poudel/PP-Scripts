standardise_data <-
function(x){
  mean_d=mean(x)
  var_d=var(x)
  
  std_data=(x - mean_d)/sqrt(var_d)
  return(std_data)
  
}
