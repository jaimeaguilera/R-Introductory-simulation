semilla1=function(a=1583458089, m=2^31-1, b=0,vecN=c(50,100,500,1000,2000,5000,10000)) {
  #Dado un cierto método y un vector N, encuentra una semilla x0 que pase
  #TODOS los tests para TODOS los números dados en vecN
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
