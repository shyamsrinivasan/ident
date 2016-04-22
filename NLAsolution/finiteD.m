% function to perform finite differenc e approximation of jacobians
function fh = finiteD
% function will only fill nonemoty spaces in the sparce Sv and consequently
% the jacobian matrix
if nargin<3
    type = 'F';
end

switch type
    case 'F'
        fh = @FD_F;        
    case 'C'
        fh = @FD_C;
    otherwise
        fprintf('Invalid option\n')
        error('FD:InvOpt','Finite Difference was not performed');        
end

function FD_C
% nothing yet to see here

function J = FD_F(Jpatt,x,fh)
% rows function evaluation
% columns state variable vector
nr = size(Jpatt,1);
% evaluate non-sparse indices
[rw,cl] = find(Jpatt);
% combine and sort on a row basis
rws = sortrows([rw cl],1);
% work on rws
Jsp = zeros(size(rws,1),1);
for irow = 1:size(Jpatt,1)
    fm = fh(irow,x);
    cls = rws(rws(:,1)==irow,2);
    for jcls = 1:length(cls)
        xmod = x;
        xmod(cls(jcls)) = xmod(cls(jcl))+eps;
        fj = fh(irow,xmod);
        Jsp(cnt) = (fj-fm)/eps;
        cnt = cnt + 1;
    end
end
J = sparse(rws(1),rws(2),Jsp,nr,nr);

% Jij = (f(x+e,p)-f(x,p))/e
% Jij = (fi([x1 x2 ... xj+e ... xn]) - fi([x1 x2 ... xj ... xn]))/e



