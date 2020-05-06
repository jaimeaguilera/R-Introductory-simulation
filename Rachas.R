Rachas=function(a) {
  n=length(a)
  if(a[1]<a[2]) t=1 else t=0  #1 si +, 0 si - 
  r=1 #número de rachas
  for(i in 2:n-1) {
   if(t==0 & a[i]<a[i+1])
    {r=r+1; t=1}
   else if(t==1 & a[i]>a[i+1])
    {r=r+1; t=0}
  }
  z=abs((r-1/3*(2*n-1))/(sqrt(1/90*(16*n-29))))
  return(z)
  print(paste(r, "rachas"))
}
