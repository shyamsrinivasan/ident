% get final equilibrium values for initial values perturbed from 
% unstable manifold
function [allxeq,ssid,allfeq] = getPequilibrium(ivals,model,pvec,options,opts,tspanf,hfig,id)

[npts,nvar] = size(ivals);

posp = ivals + repmat(1e-3*[1 1 1],npts,1);
negp = ivals - repmat(1e-3*[1 1 1],npts,1);

allxeq = zeros(nvar,npts);
allfeq = zeros(5,npts);
ha = [];
ssid = zeros(1,npts);
colorSpec = chooseColors(2,{'Green','Purple'});

for jval = 1:npts    
%     [xae1,feq1] = solveAEonly(1,posp(jval,:)',model,pvec,options);
%     [xae2,feq2] = solveAEonly(1,negp(jval,:)',model,pvec,options);      
    [~,xde1,~,feq1] = solveODEonly(1,posp(jval,:)',model,pvec,opts,tspanf);
%     [~,xde2,~,feq2] = solveODEonly(1,negp(jval,:)',model,pvec,opts,tspanf); 
    
    allxeq(1:nvar,jval) = xde1;
%     allxeq(nvar+1:end,jval) = xde2;
    allfeq(1:5,jval) = feq1;
%     allfeq(5+1:end,jval) = feq2;
    
    % equilibrium point plots    
    if xde1(1)<xde1(2) % ss2
        ssid(1,jval) = 2;
        Point.MarkerFaceColor = colorSpec{1};
        Point.MarkerEdgeColor = colorSpec{1};
%         [hfig,ha] = FIGmssEqIvalPerturbations(posp(jval,:)',xde1,2,id,hfig,ha,Point);        
    elseif xde1(1)>xde1(2) % ss1
        ssid(1,jval) = 1;
        Point.MarkerFaceColor = colorSpec{2};
        Point.MarkerEdgeColor = colorSpec{2};
%         [hfig,ha] = FIGmssEqIvalPerturbations(posp(jval,:)',xde1,2,id,hfig,ha,Point);  
    end
%     if xde2(1)<xde2(2) % ss2
%         ssid(2,jval) = 2;
%         Point.MarkerFaceColor = colorSpec{1};
%         Point.MarkerEdgeColor = colorSpec{1};
%         [hfig,ha] = FIGmssEqIvalPerturbations(negp(jval,:)',xde2,2,id,hfig,ha,Point);
%     elseif xde2(1)>xde2(2) % ss1
%         ssid(2,jval) = 1;
%         Point.MarkerFaceColor = colorSpec{2};
%         Point.MarkerEdgeColor = colorSpec{2};
%         [hfig,ha] = FIGmssEqIvalPerturbations(negp(jval,:)',xde2,2,id,hfig,ha,Point);
%     end
    
end

% collect all steady states
ss1ival = posp(ssid==1,:);
ss1eq = allxeq(:,ssid==1);
Point1.MarkerFaceColor = colorSpec{2};
Point1.MarkerEdgeColor = colorSpec{2};
ss2ival = posp(ssid==2,:);
ss2eq = allxeq(:,ssid==2);
Point2.MarkerFaceColor = colorSpec{1};
Point2.MarkerEdgeColor = colorSpec{1};

ha = [];
for iss1 = 1:size(ss1ival,1)
    [hfig,ha] =...
    FIGmssEqIvalPerturbations(ss1ival(iss1,:)',ss1eq(:,iss1),2,id,hfig,ha,Point1);  
    drawnow
end

for iss2 = 1:size(ss2ival,1)
    [hfig,ha] =...
    FIGmssEqIvalPerturbations(ss2ival(iss2,:)',ss2eq(:,iss2),2,id,hfig,ha,Point2);  
    drawnow
end





