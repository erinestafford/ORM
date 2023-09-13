function out = rf_sa(l0,I0,beta,gamma,r,c,sigma,pop,amin,amax,T)
out = get_Ttime(l0,I0,beta,gamma,r,c,sigma,pop,amin,amax)-T;
end