function [t1,t2]= get_switches(I0,beta,gamma,r,c,s,pop,T, amin, amax)
[t,~,~,a]=opSIS_Icost(I0,beta,gamma,r,c,s,pop,T, amin, amax);
t1=0;
t2=0;
temp=t(a<0.5);
if ~isempty(temp)>0
    t1=temp(1);
    temp2=t(t>temp(end));
    if ~isempty(temp2)>0
        t2=temp2(1);
    end
end
end