function [LPmax,LPmin] = solvemetLP(newmodel,prxnid)
if nargin<2
    prxnid = 0;
end

%temporary call
% newmodel = setupMetLP(model);
nmet = size(newmodel.S,1);

A = newmodel.A;
b = newmodel.b;
lb = newmodel.lb;
ub = newmodel.ub;

if ~isfield(newmodel,'cprod') && prxnid
    cprod = sparse(1,prxnid,1,1,size(A,2));
    maxMin = 1;
elseif ~isfield(newmodel,'cprod') && ~prxnid
    cprod = sparse(1,size(A,2));
    maxMin = 0;
elseif isfield(newmodel,'cprod') && ~prxnid
    if length(newmodel.cprod)==size(A,2)
        cprod = newmodel.cprod;   
    else
        n_slack = size(A,2)-length(newmodel.cprod);
        cprod = [newmodel.cprod sparse(1,n_slack)];
    end
    maxMin = 1;
end

%maximization
[x,xobj,flag,output,lambda] = cplexlp(-cprod(:),[],[],A,b,lb,ub);
if flag>0    
    LPmax.x = x;
    LPmax.obj = -xobj;
end
LPmax.flag = flag;

%minimization
if maxMin
    [x,xobj,flag,output,lambda] = cplexlp(cprod(:),[],[],A,b,lb,ub);
    if flag>0    
        LPmin.x = x;
        LPmin.obj = xobj;
    end
    LPmin.flag = flag;
else
    LPmin = struct([]);
end




    
