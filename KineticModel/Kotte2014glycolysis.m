function out = Kotte2014glycolysis
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

function dM = fun_eval(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne)
dM = zeros(4,1);
                     
%J(E, acetate)
flux(1) = kEcat.*kmrgd(1).*kmrgd(4)./(kmrgd(4)+KEacetate);
%vFbp(PEP,FBP)
ratio = 1+kmrgd(3)/KFbpFBP;
flux(3) = vFbpmax.*(ratio-1).*(ratio).^4/(ratio.^4+Lfbp*(1+kmrgd(2)./KFbpPEP));
%vEX(PEP)
flux(2) = vEXmax.*kmrgd(2)./(kmrgd(2)+KEXPEP);
%enzyme production fluxes
%E(FBP) for J (%FBP ---| Cra and Cra ---> E)
flux(4) = vemax.*(1-1./(1+(KeFBP./kmrgd(3)).^ne));

%differential equations
%enzymes
%E
dM(1) = flux(4) - d*kmrgd(1);
%PEP
dM(2) = flux(1) - flux(2);
%FBP
dM(3) = flux(2) - flux(3);
%acetate
dM(4) = 0;

function [tspan,y0,options] = init
handles = feval(Kotte2014glycolysis);

% obtain initial steady states
M = zeros(4,1);
M(1)  = 1;      % E
M(2)  = 0.01;   % PEP
M(3)  = 0.03;   % FBP
M(4) = 2;       % a.u acetate

[~,yout] = ode45(@Kotte_glycolysis,0:0.1:30,M);
y0 = yout(end,:);
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

function jac = jacobian(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne)

function jacp = jacobianp(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne)      
                     
function jacp = hessians(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne)           
                     
function jacp = hessiansp(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne)      

function jacp = der3(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne)    
                     
function jacp = der4(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne)    
                     
function jacp = der5(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne)                         




