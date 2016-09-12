function Nxeq = eqperturbation(xeq,npts,type)
if nargin<3
    type = 'rnd';
end
if nargin<2
    npts = 1;
end
nvar = size(xeq,1);

switch type
    case 'rnd'
        Nxeq = repmat(xeq,1,npts)+...
               randi([-1 1],nvar,npts).*repmat(xeq,1,npts).*...
               (2*random(makedist('Uniform'),length(xeq),npts));
        Nxeq = max(0.00001,Nxeq);
%         for i = 1:npts
%             Nxeq(:,i) = max(0.00001,Nxeq(:,1));
%         end
    case 'pos'
        Nxeq = repmat(xeq,1,npts)+...
               repmat(xeq,1,npts).*random(makedist('Uniform'),length(xeq),npts);
    case 'neg'
        Nxeq = repmat(xeq,1,npts)-...
               repmat(xeq,1,npts).*random(makedist('Uniform'),length(xeq),npts);
        Nxeq = max(0.00001,Nxeq);
%         for i = 1:npts
%             Nxeq(:,i) = max(0.00001,Nxeq(:,1));
%         end
    otherwise
        fprintf('Invalid option\nReturning replicated input\n');
        Nxeq = repmat(xeq,1,npts);
        return
end