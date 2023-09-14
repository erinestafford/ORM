function roots = get_roots(I0,beta,gamma,r,c,sigma,pop,T,amin,amax)

l0 = -1.0;
dy = 0.01;
HT = 200;
%%
lvals = [];
tts = [];
y2 = l0;
count=1;
while y2>-HT
    tt = get_Ttime(y2,I0,beta,gamma,r,c,sigma,pop,amin,amax);
    if tt<(1e6/2)
        lvals(count) = y2;
        tts(count) = tt;
        count=count+1;
    else
        break
    end
    y2=y2-dy;
    
end
%%
roots = [];
f = @(l0) rf_sa(l0,I0,beta,gamma,r,c,sigma,pop,amin,amax,T);
for i = 2:length(tts)
    if ((tts(i)<T)&&(T<tts(i-1))) ||  ((tts(i)>T)&&(T>tts(i-1)))
        lval=fzero(f,[lvals(i),lvals(i-1)]);
        roots=[roots;lval];
    end
        
end



end