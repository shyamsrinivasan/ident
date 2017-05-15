function SAxeq = getEqfromLPs(model,pvec,opts,s1,x1,contvarid)

ac = find(strcmpi(model.mets,'ac[e]'));

% get saddle node
eps1 = 1e-4;
saddle = [];
while isempty(saddle)
    [saddle,saddlepar] = getsaddlenode(s1,x1,eps1);
    eps1 = eps1*10;
end

% get eigen value and eigne vector at saddle
pvec(contvarid) = saddlepar;
model.PM(ac-length(saddle)) = saddlepar;
[~,lambda,w] = getKotteJacobian(saddle,pvec,model);

 % get 2 steady states by perturbing initial states from
                % saddle
SAxeq = [];
iter = 1;
if ~isempty(saddle)
    while iter <= 2
        delx = [1 -1;1 -1;1 -1];
        ival = saddle+1e-2*delx(:,iter);
        [~,xeq] = solveODEonly(1,ival,model,pvec,opts,1:2000);
        if ~isempty(SAxeq)
            if ~any(abs(SAxeq-repmat(xeq,1,size(SAxeq,2)))<=1e-8)
                SAxeq = [SAxeq xeq];
            end
        else
            SAxeq = xeq;
        end   
        iter = iter+1;
    end
end