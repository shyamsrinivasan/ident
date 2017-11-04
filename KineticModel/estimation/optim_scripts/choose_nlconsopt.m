function optsol = choose_nlconsopt(prob,x0,optimdata,solveropt)
if nargin<4
%     solveropt = struct('solver','ipopt','multi',1,'multi_pts',[2 2]);
    solveropt = struct('solver','ipopt','multi',0);
end

nval = size(x0,2);
if nval>1
    % use multi start search    
    optsol = multistart_search(prob,x0,optimdata,solveropt);
else
    % use single start search with a possibility for multistart
    optsol = nlconsopt(prob,x0,solveropt,optimdata);
end