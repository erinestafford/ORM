function out = bisect_ttarg(fn,a,b,tol,TTARG)
it = 1;
mid = (a+b)/2;
fa = fn(a,TTARG);
fb = fn(b,TTARG);
fmid = fn(mid,TTARG);
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
    fmid = fn(mid,TTARG);
end
out = mid;
end
