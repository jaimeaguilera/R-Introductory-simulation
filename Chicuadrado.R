Chicuadrado=function(a){
  #Proporciona el valor del estadístico chi cuadrado 
  #para un vector a de valores entre (0 y 1)
  
  #Dividimos en k=10 clases, según el primer decimal sea 0,1,..
  n=length(a)
  a=floor(10*a)  #Enteros de 0 a 10
  O=as.vector(table(factor(a,levels=0:9)))  #Vector de frecuencias
  e=n/10
  U=1/e*sum((O-e)^2)
  return(U)
}