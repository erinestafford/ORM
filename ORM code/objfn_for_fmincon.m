function objfn_val = objfn_for_fmincon(I0,amin,amax,r,sigma,beta,pop,c,gamma,Tend,t1)
t1 = abs(t1);
f = @(t,y)odeI_givent1(t,y,amin,amax,r,sigma,beta,pop,c,gamma,t1);
t = 0:.01:Tend;
[T,Y] = ode45(f,t,I0,odeset('RelTol',1e-10));
a = zeros(length(T),1);
count = 1;
for t=T'
    if t<t1
        a(count) = amin;
    else
        a(count) = amax;
    end
    count = count+1;
end
objfn_val =-trapz(T,r*a.*(pop-sigma*Y) - c*Y);
end