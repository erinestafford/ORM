function [term_times,l_vals] = find_term_times(N,target_time,beta,r,c,gamma,amax,amin,I0,pop,sigma)
    tend=1000;
    term_times = [];
    lvals=linspace(-20,-1,N);
    t=0:0.1:tend;
    
    Opt1    = odeset('Events', @(t,y)lzero(t, y));
    f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
    
    for l0 = lvals
        Y0 = [l0;I0];
        [T,~] = ode45(f,t,Y0,Opt1);
        term_times=[term_times;T(end)];
    end
    
    [~,i]=mink(abs(term_times - target_time),10);
    term_times= term_times(i);
    l_vals= lvals(i);
end
