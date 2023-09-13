function [lvals,tts] = get_term_times(l0,sigma,HT,dl)
    lvals=[];    
    tts=[];
    c=1;
    y2 = l0;
    while y2>-HT
        y2=y2-dl;
        tt = ttime_sigma(y2,sigma);
        if tt<(1e6/2)
            lvals(c) = y2;
            tts(c) = tt;
            c=c+1;
        else
            break
        end
    end
end