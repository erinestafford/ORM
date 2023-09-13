pop=1;
beta = 0.06;%0.0232/pop;%0.06;
r = 0.1;%0.39;
c = 0.2;%6.83;
gamma = 1.0/110.0;
amax = 0.9;
amin = 0.2;
I0 = 0.99;
sigma = 0.5;
T = 150;
I_term = 0.8;
%%
test=[];
lvals=-100:0.01:100;
for l=-100:0.01:100
    test=[test;V_all(I0,beta,gamma,r,c,sigma,pop,T, amin, amax, l)];
end
pairs=[];
for i=1:(length(test)-1)
    if (test(i)> I_term && test(i+1)< I_term) ||(test(i)< I_term && test(i+1)> I_term)
        pairs=[pairs;lvals(i),lvals(i+1)];
    end
end
%%
s=size(pairs);
s=s(1);
%%
for i = 1:s
    lt = secantcode(pairs(i,1),pairs(i,2),1e-4,I_term,I0,beta,gamma,r,c,sigma,pop,T, amin, amax);

    [t,I,curlambda,a]=full_eval_fbs(I0,beta,gamma,r,c,sigma,pop,T, amin, amax, lt);
    plot(t,I,linewidth=2)
    hold on
    plot(t,a,linewidth=2)
end

%%

function [I] = opSIS_fbs(I0,beta,gamma,r,c,s,pop,T, amin, amax, l_term)
test = -1;
delta = 0.001;
N = 1000;
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
curlambda(end) = l_term;
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

function [t,I,curlambda,a] = full_eval_fbs(I0,beta,gamma,r,c,s,pop,T, amin, amax, l_term)
test = -1;
delta = 0.001;
N = 1000;
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
curlambda(end) = l_term;
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

function out = V_all(I0,beta,gamma,r,c,s,pop,T, amin, amax, l_term)
I = opSIS_fbs(I0,beta,gamma,r,c,s,pop,T, amin, amax, l_term);
out = I(end);
end

function y = secantcode(a,b,epsilon,I_term,I0,beta,gamma,r,c,s,pop,T, amin, amax)
flag = -1;
V = @(l_term)I_term-V_all(I0,beta,gamma,r,c,s,pop,T, amin, amax, l_term);
Va = V(a);
Vb = V(b);
iter=0;
while(flag<0)
    if(abs(Va)>abs(Vb))
        k=a;
        a=b;
        b=k;
        k=Va;
        Va=Vb;
        Vb=k;
    end
    d = Va*(b-a)/(Vb - Va);
    b = a;
    Vb = Va;
    a = a-d;
    Va = V(a);
    if(abs(Va)<epsilon || iter>1e6)
        flag=1;
    end
    iter = iter+1;
end
y = a;
end
