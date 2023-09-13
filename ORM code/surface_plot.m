clear all; close all;

%% General surface
I_vals = linspace(0.01,0.99,100);
sigma_vals = linspace(0,1,100);
surface = zeros(100,100);
ci = 1;
for i = I_vals
    cj = 1;
    for j=sigma_vals
        surface(ci,cj) = -(1 - sigma_vals(cj)*I_vals(ci))/((1 - I_vals(ci))*I_vals(ci));
        cj=cj+1;
    end
    ci=ci+1;
end
figure()
surf(I_vals,sigma_vals,surface)
xlabel('I')
ylabel('\sigma')
zlabel('switching condition')

%% change in switching condition over one simulation

[t,I,l,a]=opSIS_Icost(0.5,0.06,1/110,0.1,0.2,0.5,1,90, 0.2, 0.9);

switch_cond = l+0.1*(1 - 0.5*I)./(0.06*(1 - I).*I);

figure()
plot(t,switch_cond, 'linewidth',3)
hold on
plot(t,zeros(1001,1), 'linewidth',3, 'linestyle','--')
xlabel('t')
ylabel('\phi')

%% Change in switching condition over time with changing sigma

H=figure();
filename = 'Switch_plot.gif';
I0_vals = 0.1:0.01:0.99999;
s_vals = linspace(0,1,100);
n=0;
for i = I0_vals
    switch_cond_var = zeros(100,1001);
    c = 1;
    for s = s_vals
        [t,I,l,a]=opSIS_Icost(i,0.06,1/110,0.1,0.2,s,1,90, 0.2, 0.9);
        switch_cond_var(c,:) = l+0.1*(1 - s*I)./(0.06*(1 - I).*I);
        c = c+1;
    end
    h =surf(t,s_vals,switch_cond_var);
    hold on
    surf(t,s_vals,zeros(100,1001))
    set(h,'LineStyle','none')
    xlabel('t')
    ylabel('\sigma')
    zlabel('\phi')
    view(-159,46)
    title('I_0 = ', i)
    colormap jet;
    pause(0.001)
    Frame = getframe(H);
    im = frame2im(Frame);
    [imind,CM] = rgb2ind(im,256);
      if n == 0
          imwrite(imind,CM,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,CM,filename,'gif','WriteMode','append');
      end
     n=n+1;
     hold off
end


