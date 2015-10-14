function newmodel = setupSlackVariables(model)

newmodel = model;

%setup slack problem
nconstr = length((newmodel.A(:,1)));
%add slack variables to all constraints
A_slack = sparse(1:nconstr,1:nconstr,1,nconstr,nconstr);

%adjust sign of slack variables to be >= 0
A_slack = repmat(sign(newmodel.Vss),1,nconstr).*A_slack;

newmodel.A = [newmodel.A A_slack];

lb_slack = zeros(size(newmodel.A,1),1);
% lb_slack(lb_slack==0) = -Inf;
ub_slack = zeros(size(newmodel.A,1),1);
ub_slack(ub_slack==0) = Inf;

newmodel.lb = [newmodel.lb;lb_slack];
newmodel.ub = [newmodel.ub;ub_slack];