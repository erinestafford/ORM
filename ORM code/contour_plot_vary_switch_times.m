%close all; clear all;

pop=100;
beta = 0.0232/pop;%0.1/pop;%0.06/pop;
r = 0.1;%0.12;%0.2;
c =0.23; %0.2;%1.41;%155/36;%6.39;
gamma = 1.0/110.0;%1.0/250;%1.0/36.0;
amax = 0.9;%0.75;
amin = 0.2;
I0 = 0.99*pop;%0.63*pop;
sigma = 0.1;%0.05;


%% vary t1 and t2

Final_time = 305;%115;
trange=0:.01:Final_time;
tplot1=0:1:Final_time;
tplot2=0:1:Final_time;
n = length(0:1:Final_time);
objf= NaN(n,n);
dM = amax*pop*beta - gamma;
dm = amin*pop*beta - gamma;

i=1;
for t1 = 0:1:Final_time
    j = 1;
    for t2 = 0:1:Final_time
        if t1<=t2
            a_vobj = zeros(length(trange),1);
            I_vobj = zeros(length(trange),1);
    
            a_vobj(trange<=t1)=amax*ones(length(a_vobj(trange<=t1)),1);
            a_vobj(trange>t1 & trange<=t2)=amin*ones(length(a_vobj(trange>t1 & trange<=t2)),1);
            a_vobj(trange>t2)=amax*ones(length(a_vobj(trange>t2)),1);
    
            I_1= @(t) dM*I0./(amax*beta*I0 + (dM - amax*beta*I0)*exp(-dM*t));
            I_2= @(t) dm*I_1(t1)./(amin*beta*I_1(t1) + (dm - amin*beta*I_1(t1))*exp(-dm*(t-t1)));
            I_3= @(t) dM*I_2(t2)./(amax*beta*I_2(t2) + (dM - amax*beta*I_2(t2))*exp(-dM*(t-t2)));


            I_vobj(trange<=t1) = I_1(trange(trange<=t1));
            I_vobj(trange>t1 & trange<=t2) = I_2(trange(trange>t1&trange<=t2));
            I_vobj(trange>t2) = I_3(trange(trange>t2));
                              
            objf(i,j) = trapz(trange,r*a_vobj.*(pop-sigma*I_vobj) - c*I_vobj);
        else
            objf(i,j) = NaN;
        end
        j = j+1;
    end
    i= i+1;
end
%% find max for plotting
[m,i]=max(transpose(objf(:)));
[i1,i2]=ind2sub(size(transpose(objf)),i);

%% Plot surface and contour
figure()
surfc(tplot1,tplot2, transpose(objf), transpose(objf), 'EdgeColor','none', 'LineWidth',1.5);
zlabel('Payoff', 'FontSize',19.8)
zlim([-2400,-800])
xlabel('t_1', 'FontSize',19.8)
ylabel('t_2', 'FontSize',19.8)
colorbar('Limits', [-2400,-800], 'FontSize',14)
ax = gca;
ax.FontSize = 18; 
% set(h,'LineStyle','none');
hold on
%add maxes
scatter3(tplot1(i1),tplot2(i2),m,50,'r','filled')
scatter3(tplot1(i1),tplot2(i2),-2400,50,'r','filled')
