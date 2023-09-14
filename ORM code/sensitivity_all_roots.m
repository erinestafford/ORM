function [out,n_switches] = sensitivity_all_roots(I0,beta,gamma,r,c,sigma,pop,T,amin,amax)
out = NaN(3,1);
n_switches=NaN(3,1);
%% find roots
roots = get_roots(I0,beta,gamma,r,c,sigma,pop,T,amin,amax);
%% run model for each root
    % get number of switches for each
    % calculate objective function value
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
options = odeset('Events', @(t,y)lzero(t, y),'RelTol',1e-8,'AbsTol',1e-10);
for l0 = 1:length(roots)
    %try stiff solver
    [t,Y] = ode23s(f,[0, T],[roots(l0),I0],options);
    [a,n_switch] = get_a(Y(:,2),Y(:,1),pop,r,sigma,beta,amin,amax);
    out(l0) = trapz(t,r*a.*(pop-sigma*Y(:,2)) - c*Y(:,2));
    n_switches(l0) = n_switch;
end
end