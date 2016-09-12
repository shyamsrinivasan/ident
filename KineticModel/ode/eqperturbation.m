function Nxeq = eqperturbation(xin,idx,npts,type)
if nargin<4
    type = 'rnd';
end
if nargin<3
    npts = 1;
end
if nargin<2
    % index of variables to be perturbed
    idx = [];
end

if ~isempty(idx)
    npvar = length(idx);    
    pval = getperturbationpoints(xin(idx),npvar,npts,type);
    Nxeq = repmat(xin,1,npts);
    Nxeq(idx,:) = pval(1:npvar,:);
else
    Nxeq = getperturbationpoints(xin,length(xin),npts,type);    
end

function pval = getperturbationpoints(xeq,npvar,npts,type)
switch type
    case 'rnd'
        pval = repmat(xeq,1,npts)+...
               randi([-1 1],npvar,npts).*repmat(xeq,1,npts).*...
               (random(makedist('Uniform'),npvar,npts));
    case 'pos'
        pval = repmat(xeq,1,npts)+...
               repmat(xeq,1,npts).*random(makedist('Uniform'),npvar,npts);
    case 'neg'
        pval = repmat(xeq,1,npts)-...
               repmat(xeq,1,npts).*random(makedist('Uniform'),npvar,npts);
    otherwise
        fprintf('Invalid option\nReturning replicated input\n');
        pval = repmat(xeq,1,npts);
        return
end
