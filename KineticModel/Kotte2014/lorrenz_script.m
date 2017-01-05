% lorrenz attractor script
% continuation on lorrenz attractor from stable equilibrium point (0,0,0)
global sys
sys.gui.pausespecial=0;  %Pause at special points 
sys.gui.pausenever=1;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

sigma = 10;
rho = 28;
beta = 8/3;
p = [sigma;rho;beta];
ap = 1;
xeq = [0;0;0];
% continuation from initial equilibrium - initialization
[x0,v0] = init_EP_EP(@KotteMATCONT,xeq,p,ap);

% MATCONT options
opt = contset;
opt = contset(opt,'VarTolerance',1e-3);
opt = contset(opt,'VarTolerance',1e-3);
opt = contset(opt,'FunTolerance',1e-3);
opt = contset(opt,'MaxNumPoints',bfpts);
opt = contset(opt,'MaxStepsize',.01);
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Eigenvalues',1);

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

if ~isempty(s1)
    data.s1 = s1;
    data.x1 = x1;
    data.f1 = f1;
    data.flux = flux1;
else
    data = [];
end