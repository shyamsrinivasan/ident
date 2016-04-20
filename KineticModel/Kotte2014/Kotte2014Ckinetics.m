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

function dM = fun_eval(t,kmrgd,model,pvec)
% pvec = [kEcat,KEacetate,...
%         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
%         vEXmax,KEXPEP,...
%         vemax,KeFBP,ne,acetate,d];
dM = zeros(length(kmrgd),1);
tfenz = find(strcmpi(model.mets,'enz[c]'));

allmc = [kmrgd;model.PM];
d = pvec.d;

% substitute with Convenience Kinetics
flux = Kotte2014Ckinetics_flux(allmc,model,pvec); 
                     
% differential equations
tfm = cellfun(@(x)strcmpi(model.mets,x),{'fdp[c]','pep[c]'},'UniformOutput',false);
tfm = cell2mat(cellfun(@(x)find(x),tfm,'UniformOutput',false));
dM(tfm) = model.S(tfm,:)*flux;

% enzymes
% E
dM(tfenz) = flux(strcmpi(model.rxns,'ENZC')) - d*allmc(tfenz);

function [tspan,y0,options] = init
handles = feval(Kotte2014Ckinetics);

% obtain initial steady states
M = zeros(3,1);
M(1)  = 1;      % E
M(2)  = 0.001;   % PEP
M(3)  = 10;   % FBP

% rhsfn = handles{2};
% rhsodefn = @(t,x)rhsfn(t,x,)

% substitute this with SUNDIALS
[~,yout] = ode45(@fun_eval,0:0.1:30,M);
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
