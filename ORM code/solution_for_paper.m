%% Uses FBS to solve bovine mastitis problem
figure()
pop_test=100;
beta_test = 0.0232/pop_test;
r_test = 0.2;%0.20;
c_test = 1.41;%155.0/110.0;
gamma_test = 1.0/110.0;
amax_test = 0.9;
amin_test = 0.2;
I0_test = 0.63*pop_test;%0.63*pop;
sigma_test = 0.05;%0.05;
[t,I,curlambda,a] = opSIS_Icost(I0_test,beta_test,gamma_test,r_test,c_test,sigma_test,pop_test,305, amin_test, amax_test);
plot(t,a, 'linewidth',3)
hold on
plot(t,I/pop_test, 'linewidth',3)
trapz(t,r*a.*(pop-sigma*I) - c*I)
        