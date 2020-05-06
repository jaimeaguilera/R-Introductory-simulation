#Estudiamos el método multiplicativo dado por a=1583458089, m=2^31-1,
#sacado de PIERRE L'ECUYER-TABLES OF LINEAR CONGRUENTIAL GENERATORS OF DIFFERENT SIZES AND GOOD LATTICE STRUCTURE

#Importamos primero las funciones necesarias (congruencial y los tests)

#0. Congruencial
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


#Para estudiar un método, creamos la funcion testcongruencial
#Dado un método y un vector N, el algoritmo devuelve un vector
#(p,p1,p2,p3) con p la proporción de semillas que
#pasan los tests para todo N, y p1,p2,p3 la proporción por separado para CADA N y CADA test (Kolmogorov, chi2 y rachas, repsectivamente).
#Las semillas son cogidas aleatoriamente, cogiendo un total de "numsemillas", que por defecto es 20.


testcongruencial=function(a=1583458089, m=2^31-1, vecN=c(50,100,500,1000,2000,5000), numsemillas=20,b=0){
  
  vecN=sort(vecN)
  vecx0=floor(runif(numsemillas)*m)  #x0 aleatorios
  Dcont=0; Ucont=0; zcont=0; cont=0
  for(i in 1:length(vecx0)){
    v=1/m*congruencial(x0=vecx0[i],a=a,m=m,b=b,N=vecN[1])
    D=Kolmogorov(v); U=Chicuadrado(v); z=Rachas(v)
    aux=0
    if(D<1.35/sqrt(vecN[1])) {Dcont=Dcont+1; aux=aux+1}
    if(U<16.92) {Ucont=Ucont+1; aux=aux+1}
    if(z<1.96) {zcont=zcont+1; aux=aux+1}
    
    if(length(vecN)>1){
      for(j in 2:length(vecN)){
        u=1/m*congruencial(x0=v[vecN[j-1]]*m,a=a,m=m,b=b,N=vecN[j]-vecN[j-1])
        u=c(v,u)
        D=Kolmogorov(u); U=Chicuadrado(u); z=Rachas(u)
        if(D<1.35/sqrt(vecN[j])) {Dcont=Dcont+1; aux=aux+1}
        if(U<16.92) {Ucont=Ucont+1; aux=aux+1}
        if(z<1.96) {zcont=zcont+1; aux=aux+1}
        v=u
      }
    }
    if(aux==3*length(vecN)) cont=cont+1
  }
  its=length(vecx0)*length(vecN)
  p=cont/length(vecx0)
  p1=Dcont/its
  p2=Ucont/its
  p3=zcont/its
  return(c(p,p1,p2,p3))}



#Estudiamos el test para N=100
c1=testcongruencial(vecN=c(1000),numsemillas=100)
print(paste("p=",c1[1]))
print(paste("p1=",c1[2]))
print(paste("p2=",c1[3]))
print(paste("p3=",c1[4]))
#Aproximadamente el 90% de las semillas pasan los tests

#Sin embargo, si tomamos varios valores de N, e.g
#vecN=c(50,100,500,1000,2000,5000) vemos que
c2=testcongruencial()
print(paste("p=",c2[1]))
print(paste("p1=",c2[2]))
print(paste("p2=",c2[3]))
print(paste("p3=",c2[4]))
#En torno al 50% de las semillas pasan los tests para todo N!!!

#Ésto nos obliga a buscar una semilla que pase los test para muchos
#valores de N, y no solo unos pocos. Lo hacemos con la función semilla1()
#Dado un cierto método y un vector N, encuentra una semilla x0 que pase
#TODOS los tests para TODOS los números dados en vecN

semilla1=function(a=1583458089, m=2^31-1, b=0,vecN=c(50,100,500,1000,2000,5000,10000)) {
 
  vecN=sort(vecN)  
  while(0<1){   #No paramos hasta que encontremos la semilla
    
    #Primer N
    x0p=floor(runif(1)*m)   #Generamos una semilla aleatoria
    v=1/m*congruencial(x0=x0p,a=a,m=m,b=b,N=vecN[1])
    D=Kolmogorov(v); U=Chicuadrado(v); z=Rachas(v)
    if(D<(1.35/sqrt(vecN[1])) & U<16.92 & z<1.96) cont=1 else cont=0
    j=2
    #Resto de N
    while(cont==j-1 & j<=length(vecN)){  #Entra sólo si ha pasado los anteriores N
      u=1/m*congruencial(x0=v[vecN[j-1]]*m,a=a,m=m,b=b,N=vecN[j]-vecN[j-1])   #Calculamos solo lo que nos falta
      u=c(v,u)    #Lo unimos a lo que ya teníamos para obtener el vector completo
      D=Kolmogorov(u); U=Chicuadrado(u); z=Rachas(u)
      if(D<(1.35/sqrt(vecN[j])) & U<16.92 & z<1.96) cont=cont+1
      j=j+1
      v=u
    }
    if(cont==length(vecN)) return(x0p)  
  }
}



print(paste("semilla=",semilla1()))
#Por defecto, se hace para vecN=c(50,100,500,1000,2000,5000,10000)).
#Podemos tomar cualquier vector de N, e.g
vN1=floor(runif(10)*4950+50) #10 aleatorios de 50 a 5000
print(paste("semillabis=",semilla1(vecN=vN1)))


#Podríamos estudiar cualquier otro método solo con especificar los parámetros en testcongruencial y semilla1

