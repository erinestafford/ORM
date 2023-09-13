
%% two switch
pop=1;
beta = 0.06;%0.0232/pop;%0.06;
r = 0.1;%0.39;
c = 0.2;%6.83;
gamma = 1.0/110.0;
amax = 0.9;
amin = 0.2;
I0 = 0.99;
sigma = 0.5;
Tend = 115;
f = @(ts)objfn_for_fmincon2(I0,amin,amax,r,sigma,beta,pop,c,gamma,Tend,ts);

xg=ga(f,2,[],[],[],[],[0,0],[Tend,Tend]);

%%
options = optimoptions('fmincon','Display','iter','Algorithm','active-set','DiffMinChange',1e-1);

%%
[x,fval,~,~,~,grad,hessian]=fmincon(f,xg,[],[],[],[],[0,0],[Tend,Tend],[],options);
display(x(1))
display(x(1)+x(2))
eig(hessian)>0

%%
[x2,fval,~,~,~,grad,hessian]=fmincon(f,[13.8,73.6],[],[],[],[],[0,0],[Tend,Tend],[],options);
display(x2(1))
display(x2(1)+x2(2))
eig(hessian)>0


%%
[x3,fval,~,~,~,grad,hessian]=fmincon(f,[0,0],[],[],[],[],[0,0],[Tend,Tend],[],options);
display(x3(1))
display(x3(1)+x3(2))
eig(hessian)>0

%%
[x4,fval,~,~,~,grad,hessian]=fmincon(f,[41.263595812189472,20.998962190309307],[],[],[],[],[0,0],[Tend,Tend],[],options);
display(x4(1))
display(x4(1)+x4(2))
eig(hessian)>0

%%
[x5,fval,~,~,~,grad,hessian]=fmincon(f,[14.779790184082659,72.582031016119970],[],[],[],[],[0,0],[Tend,Tend],[],options);
display(x5(1))
display(x5(1)+x5(2))
eig(hessian)>0
%%
format long
out_ga = objfn_for_fmincon2(I0,amin,amax,r,sigma,beta,pop,c,gamma,Tend,x); %matlab global search
out_dsoa = objfn_for_fmincon2(I0,amin,amax,r,sigma,beta,pop,c,gamma,Tend,x2); %best from dsoa

% from roots
out_r1 = objfn_for_fmincon2(I0,amin,amax,r,sigma,beta,pop,c,gamma,Tend,x3);% no switch
out_r2 = objfn_for_fmincon2(I0,amin,amax,r,sigma,beta,pop,c,gamma,Tend,x4);% non-optimal two switch root
out_r3 = objfn_for_fmincon2(I0,amin,amax,r,sigma,beta,pop,c,gamma,Tend,x5);% best root

res=[out_ga,out_dsoa,out_r1,out_r2,out_r3];
algs=["MATLAB ga", "dsoa", "root no switch","root switch 1", "root switch 2"];
[val,idx]=min(res);
algs(idx)
%%
function objfn_val = objfn_for_fmincon2(I0,amin,amax,r,sigma,beta,pop,c,gamma,Tend,ts)
t1 = abs(ts(1));
t2 = abs(ts(2));
f = @(t,y)odeI_given_t1_t2(t,y,amin,amax,beta,pop,gamma,t1,t2);
t = 0:.01:Tend;
[T,Y] = ode45(f,t,I0,odeset('RelTol',1e-10));
a = zeros(length(T),1);
count = 1;
for t=T'
    if t<=t1
        a(count) = amax;
    elseif t<=t1+t2
        a(count) = amin;
    else
        a(count)=amax;
    end
    count = count+1;
end
objfn_val =-trapz(T,r*a.*(pop-sigma*Y) - c*Y);
end

function dydt = odeI_given_t1_t2(t,y,amin,amax,beta,pop,gamma,t1,t2)
  if t<=t1
      a = amax;
  elseif t>t1&&t<=t1+t2
      a = amin;
  else
      a = amax;
  end
 
  dydt= a * beta*(pop - y)*y - gamma*y;
end
