% function to perform finite differenc e approximation of jacobians
function finiteD(model,x,type)
% function will only fill nonemoty spaces in the sparce Sv and consequently
% the jacobian matrix
if nargin<3
    type = 'F';
end
% sparsity pattern for J
[nx,nrxn] = size(model.S);
J = sparse(0,nx);
% loop through metabolites
for ix = 1:nx
    % identify all reactions involved
    rxnid = find(model.S(ix,:)~=0);
    regid = [];
    % loop through reactions
    % identify all substrates, products and regulators
    [metid,~] = ind2sub([nx length(rxnid)],find(model.S(:,rxnid)~=0));
    [regid,~] = ind2sub([nx length(rxnid)],find(model.SI(:,rxnid)~=0));    
    Jin = sparse(1,unique([metid,regid]),1,1,nx);
    J = [J;Jin];
end
spy(J)
switch type
    case 'F'
        J = FD_F(model,x,J);
        
    case 'C'
        FD_C()
    otherwise
        fprintf('Invalid option\n')
        error('FD:InvOpt','Finite Difference was not performed');        
end

function FD_C
% nothing yet to see here

function J = FD_F(model,x,Jpatt)
for ifx = 1:size(Jpatt,1)
    Jin = Jpatt(ifx,:);
end
% Jij = (f(x+e,p)-f(x,p))/e
% Jij = (fi([x1 x2 ... xj+e ... xn]) - fi([x1 x2 ... xj ... xn]))/e



