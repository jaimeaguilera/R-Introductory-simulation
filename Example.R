#Nota. Necesitamos la función de distribución del estadístico KS: "pkolm()"
#Para ello, escribir 

#install.packages("kolmim")
#library(kolmim)



#1. Estadístico Kolmogorov Smirnov
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



#2. Estadístico Chicuadrado
Chicuadrado=function(a){

#Dividimos en k=10 clases, según el primer decimal sea 0,1,..
n=length(a)
a=floor(10*a)  #Enteros de 0 a 10
O=as.vector(table(factor(a,levels=0:9)))  #Vector de frecuencias
e=n/10
U=1/e*sum((O-e)^2)
return(U)
}



#3. Rachas
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


print("EJ. 5")
#[Ej 5]
a=c(0.36, 0.16, 0.61, 0.52, 0.17, 0.88, 0.9, 0.66, 0.04, 0.93,
    0.37, 0.21, 0.10, 0.28, 0.62, 0.68, 0.78, 0.94, 0.53, 0.46)

 #Kolmogorov
 print(paste("D=",Kolmogorov(a),"< D(N=20,alpha=0.05)=0,29408. Pasa el contraste"))

 #Chi cuadrado
 #Dividimos en clases de 0 a 0.1, 0.1 a 0.2, etc
 print(paste("U1=",Chicuadrado(a),"< Chi(k-1=9,alpha=0.05)=16.919. Pasa el contraste"))

 #Si tomamos 5 clases, de 0 a 0.2, 0.2 a 0.4 etc
 e2=4
 O2=c(4,4,3,5,4)
 U2=1/e2*sum((O2-e2)^2)
 print(paste("U2=",U2,"< Chi(k-1=4,alpha=0.05)=9.8477. Pasa el contraste"))

 #Rachas
 print(paste("Z=",round(Rachas(a),2),">1.96, de N(0.1) con alpha=0.025. No pasa el contraste"))
 print("")
 

 
 
#p-valores.
congruencialpvalue=function(x0=1243821475,a=1583458089, m=2^31-1,b=0,N=1000){
  #Dado un método congruencial, calcula un vector de N números pseudoaleatorios
  #y proporciona los pvalores correspondientes a los tests de Kolmogorov, Chi2 y Rachas
  
  u=1/m*congruencial(x0=x0,a=a,b=b,m=m,N=N)
  D=Kolmogorov(u); U=Chicuadrado(u); z=Rachas(u)
  pkol=1-pkolm(D,N)
  pchi2=1-pchisq(U,df=9)
  prach=(1-pnorm(z))*2   #Multiplicamos por 2 al considerar las dos colas
  print(paste("pKol=",round(pkol,4)))
  print(paste("pchi2=",round(pchi2,4)))
  print(paste("prach=",round(prach,4)))
  return(u)
}


#Para ver algún ejemplo, primero necesitamos la función congruencial.

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


#Ejemplos
 #a) El método multiplicativo por defecto es de
 #Pierre L'ecuyer-Tables of linear congruential generators of different sizes and good lattice structure
  print("Ejemplo a)")
  u1=congruencialpvalue()
  
 #b)Coveyou-MacPherson
  print("Ejemplo b)")
  a2=5^15; b2=1; m2=2^35; x02=floor(runif(1)*m2)
  u2=congruencialpvalue(a=a2,b=b2,x0=x02,m=m2)

