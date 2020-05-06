Kolmogorov= function(a) {
  a=sort(a)
  n=length(a)
  b1=c()
  b2=c()
  for(i in 1:n) {
    b1=c(b1,i/n)
    b2=c(b2,(i-1)/n)
  }
  D=max(c(max(b1-a),max(a-b2)))
  return(D)
}