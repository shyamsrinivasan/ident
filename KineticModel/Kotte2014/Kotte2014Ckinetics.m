function out = Kotte2014Ckinetics
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

function dM = fun_eval(t,kmrgd,pvec,acetate,d)
% pvec = [kEcat,KEacetate,...
%         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
%         vEXmax,KEXPEP,...
%         vemax,KeFBP,ne,acetate,d];
dM = zeros(3,1);

% substitute with Convenience Kinetics
flux = CKinetics(model,pvec,kmrgd,[1 2 3]);
flux = Kotte_glycolysisflux(kmrgd,pvec,flux,model);
                     
% %J(E, acetate)
% flux(1) = kEcat.*kmrgd(1).*acetate./(acetate+KEacetate);
% %vFbp(PEP,FBP)
% ratio = 1+kmrgd(3)/KFbpFBP;
% flux(3) = vFbpmax.*(ratio-1).*(ratio).^3/(ratio.^4+Lfbp*(1+kmrgd(2)./KFbpPEP).^(-4));
% %vEX(PEP)
% flux(2) = vEXmax.*kmrgd(2)./(kmrgd(2)+KEXPEP);
% %enzyme production fluxes
% %E(FBP) for J (%FBP ---| Cra and Cra ---> E)
% flux(4) = vemax.*(1-1./(1+(KeFBP./kmrgd(3)).^ne));

% differential equations
tfm = cellfun(@(x)strcmpi(model.mets,x),{'fdp[c]','pep[c]'});
tfm = [find(strcmpi(model.mets,'fdp[c]'),strcmpi(model.mets
dM([1 2]) = model.S([1 2],:)*flux;
%enzymes
%E
dM(1) = flux(4) - d*kmrgd(1);
%PEP
dM(2) = flux(1) - flux(2);
%FBP
dM(3) = flux(2) - flux(3);
%acetate
% dM(4) = 0;

function [tspan,y0,options] = init
handles = feval(Kotte2014glycolysis);

% obtain initial steady states
M = zeros(3,1);
M(1)  = 1;      % E
M(2)  = 0.001;   % PEP
M(3)  = 10;   % FBP

% substitute this with SUNDIALS
[~,yout] = ode45(@Kotte_glycolysis,0:0.1:30,M);
y0 = yout(end,:);
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),...
                 'Hessians',handles(4),'HessiansP',handles(5));
tspan = [0 10];

function jac = jacobian(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne,acetate,d)

function jacp = jacobianp(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne,acetate,d)      
                     
function jacp = hessians(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne,acetate,d)           
                     
function jacp = hessiansp(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne,acetate,d)      

function jacp = der3(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne,acetate,d)    
                     
function jacp = der4(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne,acetate,d)    
                     
function jacp = der5(t,kmrgd,kEcat,KEacetate,...
                         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
                         vEXmax,KEXPEP,...
                         vemax,KeFBP,ne,acetate,d)                         
