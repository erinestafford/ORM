function [value, isterminal, direction] = lzero(T, Y)
%checks if transversality condition is satisfied
value= Y(1)<=0;
isterminal = 1;   % Stop the integration
direction  = 0;
end