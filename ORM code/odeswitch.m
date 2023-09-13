function dydt = odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma)
  if (r*(pop - sigma*y(2))+y(1)*beta*(pop-y(2))*y(2))<0
      a = amin;
  else
      a = amax;
  end
  
  dydt = zeros(2,1);
  dydt(1) = a * (r*sigma - y(1)*beta*(pop-2*y(2)))+c+gamma*y(1);
  dydt(2) = a * beta*(pop - y(2))*y(2) - gamma*y(2);
end