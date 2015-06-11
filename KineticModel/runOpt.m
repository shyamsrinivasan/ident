function [obj,flux] = runOpt(model,flux)
nflux = model.nrxn;
nmetab = model.nmetab;
nenz = length(model.EnzName);

vlb = ones(nflux,1);
vlb(vlb == 1) = -Inf;
vub = ones(nflux,1);
vub(vub == 1) = Inf;
xlb = zeros(nmetab,1);
xub = ones(nmetab,1);
xub(xub == 1) = 10;
elb = zeros(nenz,1);
eub = ones(nenz,1);
eub(eub == 1) = 100;

%Optimization using IPOPT in OPTI TB
%Objective
biomassind = strcmp('BIOMASS',model.MetabName);
objind = (model.S(biomassind,:) > 0);
cbiomass = sparse(find(objind),1,1,nflux+nmetab+nenz,1);
obj = @(x)cbiomass'*x;

%NL Constraint Bounds
cl = zeros(nflux+nmetab+nenz,1);
cu = zeros(nflux+nmetab+nenz,1);

%Variable Bounds
lb = [vlb;xlb;elb];
ub = [vub;xub;eub];

%Initial Solution
initsol = [flux;model.MC;model.EC];

%Problem Setup
options = optiset('solver','ipopt','display','iter');
opt_obj = opti('fun',obj,'nl',@KinConstr,cl,cu,'bounds',lb,ub,'x0',initsol,'options',options);
[x,fval,exitflag,info] = solve(opt_obj);


%max vbiomass
%st Sv = 0
%   vj = ConvinienceKinetics(x,e)
%   vl <= v <= vu
%   xl <= x <= xu
%   el <= e <= eu

%Non Linear Constraint
%@KinConstr;
% variables = [flux;model.MC;model.EC];
% variables = KinConstr(variables);  

% Avu = sparse(rvind,cvind,1,nflux,nflux+nmetab+nenz);
% Avl = sparse(rvind,cvind,-1,nflux,nflux+nmetab+nenz);
% Axu = sparse(rxind,cxind+nflux,1,nmetab,nflux+nmetab+nenz);
% Axl = sparse(rxind,cxind+nflux,-1,nmetab,nflux+nmetab+nenz);
% Aeu = sparse(reind,ceind+nflux+nmetab,1,nflux,nflux+nmetab+nenz);
% Ael = sparse(reind,ceind+nflux+nmetab,-1,nflux,nflux+nmetab+nenz);
    function [J] = retJac(var)
        
        [Jsatrx,Jthermrx,Jregrx,...
     relJsatrx,relJthermrx,relJregrx] = calcJacobian(subs,pruds,act,inhib,...
                                          Ksubs,Kprod,KIact,KIinb,...
                                          s_subs,s_prod,s_act,s_inhib,...
                                          numrsat(irxn),drsatsubs,drsatprod,...
                                          vact(irxn),vinhib(irxn),...
                                          gamma(irxn));
        
    end


function AX = KinConstr(var)
%Same as Conviniece Kinetics to pass as constraints to the NLP
%Defined as a nested function

metab = var(nflux+1:nflux+nmetab);
enz = var(nflux+nmetab+1:end);

jmetab = 0;
kreg = 0;

gamma = zeros(nflux,1);
vthermo = zeros(nflux,1);
numrsat = zeros(nflux,1);
drsat = zeros(nflux,1);
vsat = zeros(nflux,1);
vact = zeros(nflux,1);
vinhib = zeros(nflux,1);
vreg = zeros(nflux,1);
newflux = zeros(nflux,1);

for irxn = 1:nflux
% Initialization
    %indices for each reaction j
    subsind = model.S(:,irxn) < 0;%substrate
    prodind = model.S(:,irxn) > 0;%product
    actind = model.SI(:,irxn) > 0;%activator
    inhibind = model.SI(:,irxn) < 0;%inhibitor
    
    nsubs = length(find(subsind));
    nprod = length(find(prodind));
    nact = length(find(actind));
    ninb = length(find(inhibind));  

    %Parameters - K & KI are vectors 
    %Stoichiometric Coefficients
    %K = [Ksubs; Kprod];
    if nsubs ~= 0
        subs = metab(subsind);%Concentrations
        s_subs = -(model.S(subsind, irxn));
        Ksubs = model.K(jmetab + 1:jmetab + nsubs);
        %elements in S are -ve for substrates        
    else
        Ksubs = [];
        s_subs = [];
        subs = [];
    end    
    if nprod ~= 0
        pruds = metab(prodind);
        s_prod = model.S(prodind, irxn);
        Kprod = model.K(jmetab + nsubs + 1:jmetab + nsubs + nprod);        
    else
        Kprod = [];
        s_prod = [];
        pruds = [];
    end
    jmetab = jmetab + nsubs + nprod;
    
    %KI = [KIact; KIinb]; 
    if nact ~= 0
        act = metab(actind);
        s_act = model.SI(actind, irxn);
        KIact = model.KI(kreg + 1:kreg + nact);        
    else
        KIact = [];
        s_act = [];
        act = [];
    end
    if ninb ~= 0
        inhib = metab(inhibind);
        %elements in SI are -ve for inhbition
        s_inhib = -(model.SI(inhibind, irxn));
        KIinb = model.KI(kreg + nact + 1:kreg + nact + ninb); 
    else
        KIinb = [];
        s_inhib = [];
        inhib = [];
    end
    kreg = kreg + nact + ninb;    

    if ~isempty(subs) && ~isempty(s_subs)
        gamma_subs = prod(subs.^s_subs);
        numrsat(irxn) = prod(subs.^s_subs);
        drsatsubs = prod((1 + subs./Ksubs).^s_subs);
    else
        gamma_subs = [];
        numrsat(irxn) = 1;
        drsatsubs = 1;
    end
    if ~isempty(pruds) && ~isempty(s_prod)
        gamma_prod = prod(pruds.^s_prod);
        drsatprod = prod((1 + pruds./Kprod).^s_prod);
    else
        gamma_prod = [];
        drsatprod = 1;
    end  
    
    % Thermodynamic Contribution
    gamma(irxn) = gamma_subs/gamma_prod;
    if ~isempty(gamma(irxn))
        if model.Keq(irxn) ~= 0%avoid divide by zero error
            vthermo(irxn) = 1 - gamma(irxn)/model.Keq(irxn);
            gamma(irxn) = gamma(irxn)/model.Keq(irxn);
        end
    else
        vthermo(irxn) = 1;        
    end
    
    % Saturation Contribution   
    drsat(irxn) = drsatsubs*drsatprod;
    vsat(irxn) = numrsat(irxn)/(drsat(irxn) - 1);
    
    % Regulatory Contribution
    if ~isempty(act) && ~isempty(s_act)
        vact(irxn) = prod(((act./KIact).^s_act)./(1 + (act./KIact).^s_act));
    else
        vact(irxn) = 1;
    end
    
    if ~isempty(inhib) && ~isempty(s_inhib)
        vinhib(irxn) = prod(1./(1 + (inhib./KIinb).^s_inhib));
    else
        vinhib(irxn) = 1;
    end
    vreg(irxn) = vact(irxn)*vinhib(irxn);
    
    % Total Net Flux    
    if ~isempty(model.Vmax(irxn))
        newflux(irxn) =...
        model.Vmax(irxn)*vthermo(irxn)*vsat(irxn)*vreg(irxn);
    else
        newflux(irxn) =...
        model.Kcat(irxn)*enz(irxn)*vthermo(irxn)*vsat(irxn)*vreg(irxn);
    end
end
% var(1:nflux) = newflux;
% var(nflux+1:end) = [metab;enz];

%Constraint coefficient matrices
AS = [model.S sparse(nmetab,nmetab+nenz)];

[rvind,cvind] = find(diag(ones(nflux,1)));
[rxind,cxind] = find(diag(ones(nmetab,1)));
[reind,ceind] = find(diag(ones(nenz,1)));

%A = [Avu;Avl;Axu;Axl;Aeu;Ael];
Avar = [sparse(rvind,cvind,1,nflux,nflux+nmetab+nenz);...
        sparse(rvind,cvind,-1,nflux,nflux+nmetab+nenz);...
        sparse(rxind,cxind+nflux,1,nmetab,nflux+nmetab+nenz);...
        sparse(rxind,cxind+nflux,-1,nmetab,nflux+nmetab+nenz);...
        sparse(reind,ceind+nflux+nmetab,1,nflux,nflux+nmetab+nenz);...
        sparse(reind,ceind+nflux+nmetab,-1,nflux,nflux+nmetab+nenz)];
    
AX = AS*[newflux;metab;enz];
end
    

end












