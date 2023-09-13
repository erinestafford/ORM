function [value, isterminal, direction] = switch2(T, Y,r,sigma,beta,pop)
value= (Y(1) >= -r*(pop- sigma*Y(2))/(beta*(pop-Y(2))*Y(2))|| Y(1)>=0);
isterminal = 1;   % Stop the integration
direction  = 0;
end