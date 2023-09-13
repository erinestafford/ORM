function dydt = odeI_givent1(t,y,amin,amax,r,sigma,beta,pop,c,gamma,t1)
  if t<t1
      a = amin;
  else
      a = amax;
  end
 
  dydt= a * beta*(pop - y)*y - gamma*y;
end