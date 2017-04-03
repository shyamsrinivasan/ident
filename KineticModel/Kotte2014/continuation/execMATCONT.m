function [data,y,p] =...
execMATCONT(funhand,fluxfunhand,xeq,pvec,ap,fluxg,model,bfpts)
if nargin<8
    bfpts = 800;
end
if nargin<7
    model = [];
end
if nargin<6
    fluxg = [];
end
    
% runMATCONT
% continuation and dynamical systems analysis using MATCONT

global sys
sys.gui.pausespecial=0;  %Pause at special points 
sys.gui.pausenever=1;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

% continuation from initial equilibrium - initialization
% ap = 12; % index for parameter to be continued on     
[x0,v0] = init_EP_EP(funhand,xeq,pvec,ap);

% MATCONT options
opt = contset;
opt = contset(opt,'VarTolerance',1e-3);
opt = contset(opt,'VarTolerance',1e-3);
opt = contset(opt,'FunTolerance',1e-3);
opt = contset(opt,'MaxNumPoints',bfpts);
opt = contset(opt,'MaxStepsize',0.1);
% opt = contset(opt,'MinStepSize',0.00001);
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
if ~isempty(fluxfunhand)
    flux = fluxfunhand(xeq,x1,p,fluxg,model,pvec,ap);
else
    flux = [];
end

if ~isempty(s1)
    data.s1 = s1;
    data.x1 = x1;
    data.f1 = f1;
    data.v1 = v1;
    data.h1 = h1;
    if ~isempty(flux)
        data.flux = flux;
    else
        data.flux = [];
    end
else
    data = [];
end