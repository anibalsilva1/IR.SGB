get_sigma <- function(phi, steps){

  v <- vector(mode = "numeric", length=length(phi))

  for(i in 1:length(phi)) {
    for(j in 1:length(steps)){
      if(phi[i] >= steps[j]){
        v[i] <- v[i]+1
      }
    }
  }
  return(v)
}
