function dydt = odefcn(t,y,a,r,sigma,beta,pop,c,gamma)
  dydt = zeros(2,1);
  dydt(1) = a * (r*sigma - y(1)*beta*(pop-2*y(2)))+c+gamma*y(1);
  dydt(2) = a * beta*(pop - y(2))*y(2) - gamma*y(2);
end