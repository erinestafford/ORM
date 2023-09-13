function v_mid = shooting_switch(target_time,v1,v2,beta,r,c,gamma,amax,amin,I0,pop,sigma)
x0 = I0;
t=0:0.01:target_time;

v_mid = (v1 + v2)/2;
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
[~,Y] = ode45(f,t, [v1;x0]);
x_a = Y(end,1);
[~,Y] = ode45(f, t, [v2;x0]);
x_b = Y(end,1);
[~,Y] = ode45(f, t, [v_mid;x0]);
x_mid =Y(end,1);

numit = 1;
while abs(x_mid) > 1e-8 && numit<1000
    if x_mid(end) == 0
        break
    elseif sign(x_mid) == sign(x_a)
        v1 = v_mid;
        [~,Y] = ode45(f,t, [v1;x0]);
        x_a = Y(end,1);
    else
        v2 = v_mid;
        [~,Y] = ode45(f, t, [v2;x0]);
        x_b = Y(end,1);
    end
    v_mid = (v1 + v2)/2;

    [~,Y] = ode45(f, t, [v_mid;x0]);
    x_mid =Y(end,1);
    numit = numit+1;

end