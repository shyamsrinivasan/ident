% get final equilibrium values for initial values perturbed from 
% unstable manifold
function getPequilibrium(ivals,model,pvec,options,opts,tspanf,hfig,id)

[npts,nvar] = size(ivals);

posp = ivals + repmat(1e-3*[1 1 1],npts,1);
negp = ivals - repmat(1e-3*[1 1 1],npts,1);

allxeq = zeros(2*nvar,npts);
allfeq = zeros(2*5,npts);
ha = [];
ssid = zeros(2,npts);
colorSpec = chooseColors(2,{'Green','Purple'});

for jval = 1:npts    
%     [xae1,feq1] = solveAEonly(1,posp(jval,:)',model,pvec,options);
%     [xae2,feq2] = solveAEonly(1,negp(jval,:)',model,pvec,options);      
    [~,xde1,~,feq1] = solveODEonly(1,posp(jval,:)',model,pvec,opts,tspanf);
    [~,xde2,~,feq2] = solveODEonly(1,negp(jval,:)',model,pvec,opts,tspanf); 
    
    allxeq(1:nvar,jval) = xde1;
    allxeq(nvar+1:end,jval) = xde2;
    allfeq(1:5,jval) = feq1;
    allfeq(5+1:end,jval) = feq2;
    
    % equilibrium point plots    
    if xde1(1)<xde1(2) % ss2
        ssid(1,jval) = 2;
        Point.MarkerFaceColor = colorSpec{1};
        Point.MarkerEdgeColor = colorSpec{1};
        [hfig,ha] = FIGmssEqIvalPerturbations(posp(jval,:)',xde1,2,id,hfig,ha,Point);        
    elseif xde1(1)>xde1(2) % ss1
        ssid(1,jval) = 1;
        Point.MarkerFaceColor = colorSpec{2};
        Point.MarkerEdgeColor = colorSpec{2};
        [hfig,ha] = FIGmssEqIvalPerturbations(posp(jval,:)',xde1,2,id,hfig,ha,Point);  
    end
    if xde2(1)<xde2(2) % ss2
        ssid(1,jval) = 2;
        Point.MarkerFaceColor = colorSpec{1};
        Point.MarkerEdgeColor = colorSpec{1};
        [hfig,ha] = FIGmssEqIvalPerturbations(negp(jval,:)',xde2,2,id,hfig,ha,Point);
    elseif xde2(1)>xde2(2) % ss1
        ssid(1,jval) = 1;
        Point.MarkerFaceColor = colorSpec{2};
        Point.MarkerEdgeColor = colorSpec{2};
        [hfig,ha] = FIGmssEqIvalPerturbations(negp(jval,:)',xde2,2,id,hfig,ha,Point);
    end
    drawnow
end
