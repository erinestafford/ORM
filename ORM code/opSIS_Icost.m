function [t,I,curlambda,a] = opSIS_Icost(I0,beta,gamma,r,c,s,pop,T, amin, amax)
%FBS for our model
test = -1;
delta = 0.0001;
N = 100000;
t = linspace(0,T,N+1);
h = T/N;
h2 = h/2;

    
I = zeros(N+1,1);
oldI = zeros(N+1,1);
I(1) = I0;
    
a = zeros(N+1,1);
a(500:end) = 0.5*ones(length(a(500:end)),1);
olda = zeros(N+1,1);
olda(500:end) = 0.5*ones(length(olda(500:end)),1);
     
curlambda = zeros(N+1,1);
oldlambda = zeros(N+1,1);
it = 0;
while(test <0 && it <100000)
    olda(:) = a(:);
    oldI(:) = I(:);
    oldlambda(:) = curlambda(:);
    for i=1:1:N
        k1 = a(i)*beta*I(i)*(pop-I(i)) - gamma*I(i);
        k2 = 0.5*(a(i)+a(i+1))*beta*(I(i) + h2*k1)*(pop-(I(i) + h2*k1)) - gamma*(I(i) + h2*k1);
        k3 = 0.5*(a(i)+a(i+1))*beta*(I(i) + h2*k2)*(pop-(I(i) + h2*k2))- gamma*(I(i) + h2*k2);
        k4 = a(i+1)*beta*(I(i) + h*k3)*(pop-(I(i) + h*k3))- gamma*(I(i) + h*k3);
        I(i+1) = I(i) + (h/6.0)*(k1 + 2*k2+2*k3+k4);
    end
    for i=0:1:N-1
        j = N+1 - i;
        k1 = a(j)*(r*s - curlambda(j)*beta*(pop-2*I(j))) + c + gamma*curlambda(j);
        k2 = 0.5*(a(j)+a(j-1))*(r*s - (curlambda(j) - h2*k1)*beta*(pop-(I(j)+I(j-1)))) + c+ gamma*(curlambda(j) - h2*k1);
        k3 =0.5*(a(j)+a(j-1))*(r*s - (curlambda(j) - h2*k2)*beta*(pop-(I(j)+I(j-1)))) + c+ gamma*(curlambda(j) - h2*k2);
        k4 = a(j-1)*(r*s - (curlambda(j) - h*k3)*beta*(pop-2*I(j-1))) + c+ gamma*(curlambda(j) - h*k3);
        curlambda(j-1)=curlambda(j)-(h/6.0)*(k1+2*k2+2*k3+k4);

    end
    anew = zeros(N+1,1);
    for i =1:1:N+1
        phi = r*(pop-s*I(i))+curlambda(i)*beta*(pop-I(i))*I(i);
        if phi <= 0
                anew(i) = amin;
        elseif phi>=0
                anew(i) = amax;
        end
    end
    a = 0.5*(a+anew);
        
    temp1=delta*sum(abs(a))-sum(abs(olda-a));
    temp2=delta*sum(abs(I))-sum(abs(oldI-I));
    temp3=delta*sum(abs(curlambda))-sum(abs(oldlambda-curlambda));

    test = min([temp1,temp2,temp3]);
    it = it+1;
end
end