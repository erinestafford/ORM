%%
pop=1;
beta = 0.06;%0.0232/pop;%0.06;
r = 0.1;%0.39;
c = 0.2;%6.83;
gamma = 1/110;
amax = 0.9;
amin = 0.2;
I0 = 0.99;

ds = 0.1;
smin=0.4;
smax=0.7;

l0 = -1.0;
HT = 20;
T=50.0:10:300;

dl = 0.1;
sigmas = linspace(smin,smax,20);
f = @(t,y)odeswitch(t,y,amin,amax,r,sigma,beta,pop,c,gamma);
all_roots=zeros(length(T),20,3);
k=1;
for TTARG=T
roots=zeros(20,3);
j = 1;
for sigma=sigmas
    [lvals,tts] = get_term_times(l0,sigma,HT,dl);
    c=1;
    roots_temp = [0,0,0];
    for i = 2:length(tts)
    if ((tts(i)<TTARG)&&(TTARG<tts(i-1))) ||  ((tts(i)>TTARG)&&(TTARG>tts(i-1)))
        lval=fzero(@(x)rf_sigma(x,TTARG,sigma),[lvals(i),lvals(i-1)]);
        roots_temp(c)=lval;
        c=c+1;
    end
    end
    
    roots(j,:)=roots_temp;
    j = j+1;
  
    
    
end

all_roots(k,:,:)=roots;
k=k+1;
end

%%

[Sval,Tvals]= meshgrid(sigmas,T);
root_plot = all_roots(:,:,1);
root_plot(root_plot==0)=nan;
surfc(Sval,Tvals,root_plot)
hold on
%contourf(Sval,Tvals,root_plot)

root_plot = all_roots(:,:,2);
%root_plot(root_plot==0)=nan;
surfc(Sval,Tvals,root_plot)
%contourf(Sval,Tvals,root_plot)
root_plot = all_roots(:,:,3);
root_plot(root_plot==0)=nan;
surfc(Sval,Tvals,root_plot)
%contourf(Sval,Tvals,root_plot)
xlabel('\sigma')
ylabel('T')
zlabel('roots')
shading interp
