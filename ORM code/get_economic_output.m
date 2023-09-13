function [out,l0] = get_economic_output(I0,beta,gamma,r,c,s,pop,T, amin, amax)
[t,I,L,a]=opSIS_Icost(I0,beta,gamma,r,c,s,pop,T, amin, amax);
out = trapz(t,r*a.*(pop-s*I) - c*I);
l0=L(1);
end