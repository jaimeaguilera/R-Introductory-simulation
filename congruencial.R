congruencial=function(x0,a,b=0,m,N=10,per=0,plt=0) {
  #método y gráfica
  z=c()
  x=x0
  x1=(a*x+b)%%m
  if(plt==1) plot(x,x1,xlim=c(0,m),ylim=c(0,m))
  x=x1
  z=c(z,x)
  for(i in 2:N) {
    x1=(a*x+b)%%m
    if(plt==1) points(x,x1)
    x=x1
    z=c(z,x)
  }
  #periodo
  if(per==1) {
   if(N<m)
    for(i in N+1:m) {
      x=(a*x+b)%%m
    }
   p=1
   xm=x
   x=(a*x+b)%%m
    while(x!=xm){
      x=(a*x+b)%%m
      p=p+1
    }
   print(paste("periodo",p))}
  
  return(z)
}