function [data,y,p] = execMATCONT(xeq,pvec,ap,fluxg,model)
% runMATCONT
% continuation and dynamical systems analysis using MATCONT

global sys
sys.gui.pausespecial=0;  %Pause at special points 
sys.gui.pausenever=1;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

% continuation from initial equilibrium - initialization
% ap = 12; % index for parameter to be continued on     
[x0,v0] = init_EP_EP(@KotteMATCONT,xeq,pvec,ap);

% MATCONT options
opt = contset;
opt = contset(opt,'VarTolerance',1e-3);
opt = contset(opt,'VarTolerance',1e-3);
opt = contset(opt,'FunTolerance',1e-3);
opt = contset(opt,'MaxNumPoints',6000);
opt = contset(opt,'MaxStepsize',.01);
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Eigenvalues',1);
% opt = contset(opt,'Backward',1);

% Equilibrium Continuation
[x1,v1,s1,h1,f1] = cont(@equilibrium,x0,v0,opt); 

% separation of variable and parameter vectors
if ~isempty(s1)    
    y = x1(1:length(xeq),:);
    p = x1(length(xeq)+1:end,:);    
else
    y = [];
    p = [];
end

% calculation  of fluxes
ac = find(strcmpi(model.mets,'ac[e]'));
flux1 = zeros(length(fluxg),size(x1,2));
if ~isempty(x1)
    for icp = 1:size(x1,2)
        pvec(ap) = p(icp);
        model.PM(ac-length(xeq)) = p(icp);
        flux1(:,icp) = Kotte_givenFlux([x1(1:length(xeq),icp);model.PM],pvec,model);
    end
end

if ~isempty(s1)
    data.s1 = s1;
    data.x1 = x1;
    data.f1 = f1;
    data.flux = flux1;
else
    data = [];
end