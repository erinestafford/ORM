data = load("data_from_cusp_attempt.mat");
roots=data.all_roots;
sigmas = data.Sval;
Ts = data.Tvals;
roots(roots==0)=nan;
figure()
scatter3(sigmas(1,:),Ts(1,:),roots(:,:,1),"bo",'filled')
hold on
scatter3(sigmas(1,:),Ts,roots(:,:,2),"ro",'filled')
scatter3(sigmas(1,:),Ts,roots(:,:,3),"go",'filled')

%%
figure()
scatter(sigmas(1,:),roots(10,:,1),"ro",'filled')
hold on
scatter(sigmas(1,:),roots(10,:,2),"ro",'filled')
scatter(sigmas(1,:),roots(10,:,3),"ro",'filled')
