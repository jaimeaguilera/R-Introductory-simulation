VonNeumann = function(x0="3434",N=10) {
  if(class(x0)=="numeric") stop("x0 no va entre comillas")
  n=nchar(x0)/2
  if (n%%1 != 0) stop("x0 no tiene 2n dígitos")
  for(i in 1:N) {
    x0=as.integer(x0)^2
    while (nchar(x0)<4*n) {x0=paste0("0",x0)}
    x0=substring(x0,n+1,3*n)
    print(x0)
  }
}	
