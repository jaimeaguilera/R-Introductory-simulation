congruencialpvalue=function(x0=1243821475,a=1583458089, m=2^31-1,b=0,N=1000){
  #Dado un método congruencial, calcula un vector de N números pseudoaleatorios
  #y proporciona los pvalores correspondientes a los tests de Kolmogorov, Chi2 y Rachas
  
  u=1/m*congruencial(x0=x0,a=a,b=b,m=m,N=N)
  D=Kolmogorov(u); U=Chicuadrado(u); z=Rachas(u)
  pkol=1-pkolm(D,N)
  pchi2=1-pchisq(U,df=9)
  prach=(1-pnorm(z))*2   #Multiplicamos por 2 al considerar las dos colas
  print(paste("pKol=",pkol))
  print(paste("pchi2=",pchi2))
  print(paste("prach=",prach))
  return(u)
}