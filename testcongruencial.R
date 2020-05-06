testcongruencial=function(a=1583458089, m=2^31-1, vecN=c(50,100,500,1000,2000,5000), numsemillas=20,b=0){
#Dado un método y un vector N, el algoritmo devuelve un vector
#(p,p1,p2,p3) con p la proporción de semillas que
#pasan los tests para todo N, y p1,p2,p3 la proporción por separado para CADA N y CADA test.
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


