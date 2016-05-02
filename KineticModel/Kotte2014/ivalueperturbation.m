function Nxeq = ivalueperturbation(lb,ub,npts)
if nargin<3
    npts = 1000;
end
if length(lb)~=length(ub)
    error('ivalpertb:sizemis','Size mismatch between lower bound and upper bound vectors\n');
end
Nxeq = repmat(lb,1,npts) + (repmat(ub,1,npts)-repmat(lb,1,npts)).*...
       random(makedist('Uniform'),length(lb),npts);
   
for i = 1:npts
    Nxeq(:,i) = max(0.00001,Nxeq(:,i));
end