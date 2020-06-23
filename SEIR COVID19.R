
#PRÁCTICA OPCIONAL.
#JAIME AGUILERA DÍAZ.
#MODELO SEIR ESTOCÁSTICO CON ALGUNAS MODIFICACIONES.
#15 DE JUNIO DE 2020.


#El objetivo es estudiar la propagación de un brote según un modelo SEIR, en el que se le ha añadido 
#un compartimento adicional de individuos fallecidos (f). Según el SEIR, dados (s,e,i,r), teníamos que
#P(s,e,i-1,r+1)=g*i,  con g la tasa de recuperaciones. En nuestro modelo, sin embargo, distinguiremos dos 
#casos, y dados (s,e,i,r,f), tendremos P(s,e,i-1,r+1,f)=g*i, P(s,e,i-1,r,f+1)=d*i, con
#d la tasa de defunciones. 

#También incluimos la posibilidad de tomar medidas, en el caso de que el brote alcance determinada
#magnitud, a fin de parar el brote. Las medidas se traducen en que disminuye la tasa de contagio b.
#En nuestro modelo, en el caso de que haya muchos infectados y pueda colapsar el sistema de salud,
#que ocurre en, digamos, i=N/10, se toma b=0.95*(g+d), de forma que la probabilidad de 
#que se produzca una nueva exposición es menor que la probabilidad de se produzca una recuperación
#o una muerte, con lo que la enfermedad tiende a remitir. 


SEIRF=function(medidas,N,a,b,g,d,s0,e0,i0,M) {
  
  #Dada una población N y una enfermedad con a,b,g,d las tasas de exposición, contactos, 
  #recuperaciones y defunciones, respectivamente, SEIRF simula la evolución del brote
  #proporcionando s,e,i,r,f, el número de susceptibles, expuestos, infectados, recuperados
  #y fallecidos, respectivamente, en cada suceso k.
  
  
  if (s0+e0+i0>N) stop("Condiciones iniciales incorrectas")
  f0=0                  #0 fallecidos inicialmente
  r0=N-s0-e0-i0         #inmunizados en el inicio del brote
  a=a/N; b=b/N; g=g/N; d=d/N   #a,b,g,d<1. Así nos aseguramos de que las probabilidades sean <1.
  s=s0; e=e0; i=i0; r=r0; f=f0
  k=0
  D=c(k,s,e,i,r,f)
  while(k<M & (e!=0 | i!=0)) {
    if(medidas=="si") {   
      if(i>N/10) {      
        b=0.95*(g+d)   #Si hay muchos infectados tomo medidas
        medidas="no"}
    } 
    U=runif(1)         #Simulo un suceso
    if (U<b*i*s/N) {s=s-1; e=e+1
    } else if (U<b*i*s/N+a*e) {e=e-1; i=i+1
    } else if (U<b*i*s/N+a*e+g*i) {i=i-1; r=r+1
    } else if (U<b*i*s/N+a*e+g*i+d*i) {i=i-1; f=f+1}
    D=c(D,k,s,e,i,r,f)
    k=k+1
  }
  dim(D)=c(6,k+1)
  if(e!=0 | i!=0) {print("El brote no ha acabado")}
  return(D)
}

#En lo que sigue, tomamos una población de N=300 individuos. Cuanto mayor es la población
#menores son los efectos estocásticos, y a N suficientemente grande, se puede usar
#un SEIR determinista en su lugar. 

#PROPAGACIÓN DE BROTES

#Si una enfermedad infecciosa tiene b mayor que g+d, es decir, es muy contagiosa y la recuperación
#o el fallecimiento es lento, la enfermedad tiende a propagarse al principio de forma exponencial
#hasta que el numero de susceptibles baja lo suficiente, frenando la propagación.
#Se tiende hacia un estado (s,0,0,f,r), con s pequeño (casi todo el mundo contrae la enfermedad).
#Notar que, en gral., el estado (s,e,0) es absorbente.

par(mfrow=c(2,2))  #Para mostrar 4 gráficas a la vez. 
#Te recomiendo poner la ventana de las gráficas bastante grande para que se vea bien

D=SEIRF(medidas="no", N=300, a=1,b=0.9,g=0.25,d=0.05,s0=295,e0=0,i0=5, M=50000)
plot(D[1,],D[2,], ylab="Susceptibles (s)", xlab="Sucesos", type="l", col="blue")
plot(D[1,],D[3,], ylab="Expuestos (e)", xlab="Sucesos", type="l", col="orange")
plot(D[1,],D[4,], ylab="Infectados (i)", xlab="Sucesos", type="l", col="red")
plot(D[1,],D[6,], ylab="Fallecidos (f)", xlab="Sucesos", type="l", col="black")


#Si una enfermedad infeccionsa tiene b menor que g+d, el brote tiende a remitir

D=SEIRF(medidas="no", N=300, a=1,b=0.25,g=0.25,d=0.05,s0=270,e0=0,i0=30,M=50000)
plot(D[1,],D[2,], ylab="Susceptibles (s)", xlab="Sucesos", type="l", col="blue")
plot(D[1,],D[3,], ylab="Expuestos (e)", xlab="Sucesos", type="l", col="orange")
plot(D[1,],D[4,], ylab="Infectados (i)", xlab="Sucesos", type="l", col="red")
plot(D[1,],D[6,], ylab="Fallecidos (f)", xlab="Sucesos", type="l", col="black")


#En la práctica, lo que hemos visto con el COVID es que, al principio, la enfermedad se propagaba con 
#fuerza, de forma exponencial, hasta que se tomaron medidas que hicieron un contagio más controlado.
#LLevamos esto a nuestro modelo.

D=SEIRF(medidas="si", N=300, a=1,b=0.9,g=0.25,d=0.05,s0=295,e0=0,i0=5, M=50000)
plot(D[1,],D[2,], ylab="Susceptibles (s)", xlab="Sucesos", type="l", col="blue")
plot(D[1,],D[3,], ylab="Expuestos (e)", xlab="Sucesos", type="l", col="orange")
plot(D[1,],D[4,], ylab="Infectados (i)", xlab="Sucesos", type="l", col="red")
plot(D[1,],D[6,], ylab="Fallecidos (f)", xlab="Sucesos", type="l", col="black")


#Si, tras ello, relajamos un poco las medidas y hay rebrote, al haber inmunidad la segunda oleada es más suave

tf=max(D[1,])
s0_=D[2,tf]
f0_=D[6,tf]
D=SEIRF(medidas="no", N=300-f0_, a=1,b=0.8,g=0.25,d=0.05,s0=s0_,e0=0,i0=5, M=50000)
plot(D[1,],D[2,], ylab="Susceptibles (s)", xlab="Sucesos", type="l", col="blue")
plot(D[1,],D[4,], ylab="Infectados (i)", xlab="Sucesos", type="l", col="red")
plot(D[1,],D[5,], ylab="Recuperados (r)", xlab="Sucesos", type="l", col="green")
plot(D[1,],D[6,], ylab="Nuevos fallecidos (f)", xlab="Sucesos", type="l", col="black")

#De hecho, cuando la inmunidad es suficientemente alta (inmunidad de rebaño),
#la enfermedad acaba desapareciendo. 



#ESTIMACIÓN DE PARÁMETROS.

#Estimamos por simulación el número máximo de infectados del brote, el número total de fallecidos,
#y sus intervalos de confianza para los casos sin medidas y con medidas. También damos un
#histograma de frecuencias para los datos simulados.
#Llamamos iM al maximo de infectados de cada brote y tf al tiempo de fin del brote. 
#f son los fallecidos de cada brote. alpha nos da un intervalo de confianza al 100*(1-alpha)%. 
#l_iM es la longitud maxima del intervalo de confianza para los máximos infectados.
#Pondremos un subíndice m para en el caso con medidas.


estimacion=function(alpha,l_iM, l_iMm, N,a,b,g,d,s0,e0,i0,M) {
  
  
  print("Brote sin medidas:")
  
  z_alpha=qnorm(1-alpha/2)
  length_iM=l_iM+1
  j=1
  iM=c()
  f=c()
  
  while(length_iM>l_iM) {
    D=SEIRF(medidas="no",N,a,b,g,d,s0,e0,i0,M)
    iM=c(iM,max(D[4,]))
    tf=max(D[1,])
    f=c(f,D[6,tf])
    if(j>10) {
      mu_iM=sum(iM)/j
      S_iM=sqrt(1/(j-1)*sum((iM-mu_iM)^2))
      length_iM=2*z_alpha*S_iM/sqrt(j)
    }
    j=j+1
  }
  #Estimación número de fallecidos al 95%
  mu_f=sum(f)/j
  S_f=sqrt(1/(j-1)*sum((f-mu_f)^2))
  length_f=2*z_alpha*S_f/sqrt(j)
  
  print(paste("iM=",round(mu_iM,1),"con intervalo de confianza al 95% de (",round(mu_iM-length_iM/2,1),",",round(mu_iM+length_iM/2,1),")"))
  print(paste("f=",round(mu_f,1),"con intervalo de confianza al 95% de (",round(mu_f-length_f/2,1),",",round(mu_f+length_f/2,1),")"))
  print(paste(j,"simulaciones"))
  hist(iM, main="", xlab="Pico de infectados (iM)", ylab="Frecuencia", col="blue")
  hist(f, main="", xlab="Fallecidos (f)", ylab="Frecuencia",col="azure3")
  
  
  #Ahora tomando medidas cuando i=N/10 (llamo iMm)
  
  print("Brote con medidas:")
  
  length_iMm=l_iMm+1
  j=1
  iMm=c()
  fm=c()
  
  while(length_iMm>l_iMm) {
    D=SEIRF(medidas="si",N,a,b,g,d,s0,e0,i0,M)
    iMm=c(iMm,max(D[4,]))
    tfm=max(D[1,])
    fm=c(fm,D[6,tfm])
    if(j>10) {
      mu_iMm=sum(iMm)/j
      S_iMm=sqrt(1/(j-1)*sum((iMm-mu_iMm)^2))
      length_iMm=2*z_alpha*S_iMm/sqrt(j)
    }
    j=j+1
  }
  mu_fm=sum(fm)/j
  S_fm=sqrt(1/(j-1)*sum((fm-mu_fm)^2))
  length_fm=2*z_alpha*S_fm/sqrt(j)
  
  
  print(paste("iMm=",round(mu_iMm,1),"con intervalo de confianza al 95% de (",round(mu_iMm-length_iMm/2,1),",",round(mu_iMm+length_iMm/2,1),")"))
  print(paste("fm=",round(mu_fm,1),"con intervalo de confianza al 95% de (",round(mu_fm-length_fm/2,1),",",round(mu_fm+length_fm/2,1),")"))
  print(paste(j,"simulaciones"))
  hist(iMm, main="", xlab="Pico de infectados(iMm)", ylab="Frecuencia", col="cyan")
  hist(fm, main="",xlab="Fallecidos (fm)", ylab="Frecuencia",col="azure4")
}



#Estimamos los parámetros para un intervalo de confianza al 95% (alpha=0.05), con l_iM=10, l_iMm=5.

estimacion(alpha=0.05,l_iM=10, l_iMm=4, N=300,a=1,b=0.9,g=0.25,d=0.05,s0=295,e0=0,i0=5,M=50000)

#OJO, suele tardar unos 30s en ejecutarse. 

#Se obtiene que el numero máximo de infectados y el número de defunciones es aproximadamente la mitad
#cuando se toman medidas. Vemos así que las restricciones encaminadas a disminuir la tasa de
#contagio b son determinantes para el control del virus, para evitar el colapso del sistema
#sanitario, y para reducir el número de defunciones.
