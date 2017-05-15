function [par] = defparval(nterms,par)
    if nargin < 2
      par = zeros(nterms,1);
    else
      par = [par;zeros(nterms,1)];
    end
    par(par == 0) = 1;
end