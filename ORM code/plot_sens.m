%% Needs to be adapted by number of switches

out_I0_no_zeros = out_I0;
out_I0_no_zeros(out_I0_no_zeros==0)=nan;

figure()
plot(I0_values,out_I0_no_zeros(:,1),'--r', 'linewidth',2)
hold on
plot(I0_values,out_I0_no_zeros(:,2),'-b', 'linewidth',2)
plot(I0_values,out_I0_no_zeros(:,3),'o-g', 'linewidth',2)
xlabel("I_0")
ylabel("Economic Output")

%%
out_sigma_no_zeros = out_sigma;
out_sigma_no_zeros(out_sigma_no_zeros==0)=nan;

figure()
plot(s_values,out_sigma_no_zeros(:,1),'--r', 'linewidth',2)
hold on
plot(s_values,out_sigma_no_zeros(:,2),'-b', 'linewidth',2)
plot(s_values,out_sigma_no_zeros(:,3),'o-g', 'linewidth',2)
xlabel("\sigma")
ylabel("Economic Output")

%%
out_T_no_zeros = out_T;
out_T_no_zeros(out_T_no_zeros==0)=nan;

figure()
plot(T_values,out_T_no_zeros(:,1),'--r', 'linewidth',2)
hold on
plot(T_values,out_T_no_zeros(:,2),'-b', 'linewidth',2)
plot(T_values,out_T_no_zeros(:,3),'o-g', 'linewidth',2)
xlabel("\sigma")
ylabel("Economic Output")