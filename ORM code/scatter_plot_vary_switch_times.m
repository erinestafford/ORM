%close all; clear all;
pop=100.0;
beta = 0.0232/pop;
r = 0.1;%0.11;%0.1 0.1;
c = 0.2;%0.2;%0.23 0.2;
gamma = 1.0/110.0;
amax = 0.9;
amin = 0.2;
sigma = 0.1;
I0 = 0.99*pop;
Final_time = 305;

%check [0,0], [128.5658,135.5019], [33.5039, 207.4193]
% pop=100;
% beta = 0.0232/pop;
% r = 0.20;
% c = 155.0/110.0;
% gamma = 1.0/110.0;
% amax = 0.9;
% amin = 0.2;
% I0 = 0.63*pop;
% sigma = 0.05;
% Final_time = 305;

%% vary t1 and t2
figure()
num_test =100000;
switch_vals = zeros(num_test,2);
switch_vals(1,:) = [0,0];
switch_vals(2,:) = [128.5501,135.4921];%[90.7823,154.5812];%[112.3208,142.7259];
switch_vals(3,:)= [33.5039,207.4193];%[64.2889,173.8194];%[45.8627, 191.8329];
for i =4:num_test
    switch_vals(i,:)=sort(randi(Final_time,1,2));
end
%Plot objective functional for general SI
objf = zeros(num_test,1);
trange=linspace(0.1,Final_time,10001);

dM = amax*pop*beta - gamma;
dm = amin*pop*beta - gamma;

count=1;
for i = 1:length(switch_vals)
    t1 = switch_vals(i,1);
    t2 = switch_vals(i,2);
    a_vobj = zeros(length(trange),1);
    I_vobj = zeros(length(trange),1);
    
    a_vobj(trange<t1)=amax*ones(length(a_vobj(trange<t1)),1);
    a_vobj(trange>=t1 & trange<t2)=amin*ones(length(a_vobj(trange>=t1 & trange<t2)),1);
    a_vobj(trange>=t2)=amax*ones(length(a_vobj(trange>=t2)),1);
    
    I_1= @(t) dM*I0./(amax*beta*I0 + (dM - amax*beta*I0)*exp(-dM*t));
    I_2= @(t) dm*I_1(t1)./(amin*beta*I_1(t1) + (dm - amin*beta*I_1(t1))*exp(-dm*(t-t1)));
    I_3= @(t) dM*I_2(t2)./(amax*beta*I_2(t2) + (dM - amax*beta*I_2(t2))*exp(-dM*(t-t2)));


    I_vobj(trange<t1) = I_1(trange(trange<t1));
    I_vobj(trange>=t1 & trange<t2) = I_2(trange(trange>=t1&trange<t2));
    I_vobj(trange>=t2) = I_3(trange(trange>=t2));
                              
    objf(count) = (Final_time/length(trange))*trapz(r*a_vobj.*(pop-sigma*I_vobj) - c*I_vobj);
    count=count+1;
end

scatter3(switch_vals(:,1),switch_vals(:,2),objf,50,objf,'filled')
hold on
% scatter3(switch_vals(1,1),switch_vals(1,2),objf(1),100,'r','filled')
% scatter3(switch_vals(2,1),switch_vals(2,2),objf(2),100,'r','filled')
% scatter3(switch_vals(3,1),switch_vals(3,2),objf(3),100,'r','filled')
[m,i]=max(objf);
scatter3(switch_vals(i,1),switch_vals(i,2),m,100,'r','filled')
zlabel('Objective Function Value')
xlabel('t1')
ylabel('t2')

colorbar()