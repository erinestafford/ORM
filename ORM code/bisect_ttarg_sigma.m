function out = bisect_ttarg_sigma(fn,a,b,tol,TTARG,sigma)
it = 1;
mid = (a+b)/2;
fa = fn(a,TTARG,sigma);
fb = fn(b,TTARG,sigma);
fmid = fn(mid,TTARG,sigma);
while (it<1e6 && abs(fmid)>tol)
    it=it+1;
    if sign(fmid)==sign(fa)
        a = mid;
        fa=fmid;
    else
        b = mid;
        fb = fmid;
    end
    mid = (a+b)/2;
    fmid = fn(mid,TTARG,sigma);
end
out = mid;
end