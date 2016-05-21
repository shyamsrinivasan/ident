% runMATCONT
% continuation and dynamical systems analysis using MATCONT
global sys
sys.gui.pausespecial=0;  %Pause at special points 
sys.gui.pausenever=1;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

% continuation from initial equilibrium - initialization
ap = 12; % index for parameter to be continued on     
[x0,v0] = init_EP_EP(@KotteMATCONT,xeq,pvec,ap);

% MATCONT options
opt = contset;
opt = contset(opt,'VarTolerance',1e-3);
opt = contset(opt,'VarTolerance',1e-3);
opt = contset(opt,'FunTolerance',1e-3);
opt = contset(opt,'MaxNumPoints',500);
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
end
