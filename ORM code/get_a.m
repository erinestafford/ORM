function [a,num_switch,switch_ind] = get_a(I,L,pop,r,sigma,beta,amin,amax)
num_switch = 0;
a = zeros(length(I),1);
switch_ind=[];
    for j = 1:length(I)
        if (r*(pop - sigma*I(j))+L(j)*beta*(pop-I(j))*I(j))<0
            a(j) = amin;
        else
            a(j) = amax;
        end
        
        if j>1 && a(j) ~= a(j-1)
            num_switch =num_switch +1;
            switch_ind = [switch_ind;j];
        end
    end
end