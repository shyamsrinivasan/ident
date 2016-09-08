function pvec = buildmodels(model,pvec,mc,rxn_add,rxn_excep)
if nargin < 5
    rxn_excep = {};
end
if nargin<4
    rxn_add = {};
end

model.rxn_add = rxn_add;
model.rxn_excep = rxn_excep;

% reactions to consider for kinetics other than Vind
Vind = addToVind(model.rxns,model.Vind,rxn_add,rxn_excep);

% metabolites that do not affect thermodynamic equilibrium  
he = find(strcmpi(model.mets,'h[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));

nrxn = length(Vind);
ntrxn = model.nt_rxn;
ntmet = model.nt_metab;

% backup known parameters from pvec
newp = struct();
newp.K = pvec.K;
newp.Kind = sparse(ntmet,ntrxn);
newp.Vmax = pvec.Vmax;
newp.kfwd = pvec.kfwd;
newp.kbkw = pvec.krev;

pvec.Kin = sparse(ntmet,ntrxn);
S = model.S;
check = zeros(ntrxn,1);
flux = zeros(ntrxn,1);

for irxn = 1:nrxn
    sbid = S(:,Vind(irxn))<0;    
    prid = S(:,Vind(irxn))>0;        
    % remove water
    sbid(h2o) = 0;
    prid(h2o) = 0;      
    % remove protons
    sbid([he hc]) = 0;
    prid([he hc]) = 0;
    
    % no parameters for cofactors - assumed abundant 
    % cofactros are assumed as compensated species
    % hence   
    rerun = 0;
    Ksbackup = pvec.K(sbid,Vind(irxn));
    Kpbackup = pvec.K(prid,Vind(irxn));
    kfwdbkup = pvec.kfwd(Vind(irxn));
    kbkwbkup = pvec.krev(Vind(irxn));
    while check(Vind(irxn))~=1 
%         if any(Ksbackup==1)||any(Kpbackup==1)
            % sampling of parameters needs to be recursive until check (see below) is 1
        pvec = estimateKm(pvec,sbid,prid,mc,Ksbackup,Kpbackup,Vind(irxn),rerun);
    
    % forward and backward catalytic rates
    % kfwd and kbkw
    % kfwd or kbkw is sampled basedon reaction directionality from FBA for
    % thermodynamic consistency
    % sampling done only for unknown values
%     fprintf('%s \t delG = %3.6g \t Vflux = %3.6g\t',model.rxns{Vind(irxn)},...
%              pvec.delGr(Vind(irxn)),...
%              model.Vss(Vind(irxn)));  
%     if ~any(strcmpi(model.rxns{Vind(irxn)},'ATPS4r'))
        % ATPsynthase does not strictly obey Briggs Haldane
        pvec = samplekcat(model,pvec,sbid,prid,Vind(irxn),mc,kfwdbkup,kbkwbkup,rerun);
%     end    
        pvec.Vmax(Vind(irxn)) = 1;
        pvec.Vmax(model.Vss==0) = 0;
    
        % check for vss and delGr direction
        check = checkthermo(@CKinetics,check);
        rerun = 1;
        % if not satisfied => resample -> use a while loop?
%         else
%             fprintf('Parameter estimation entering an infinite loop for %s\n',model.rxns(Vind(irxn)))
%         end
    end
end

% other reactions 
% transport reactions x[e] <==> x[c]
new_excep = union(rxn_excep,model.rxns(Vind));
Vex = addToVind(model.rxns,model.Vex,[],new_excep);
for irxn = 1:length(Vex)
    nmet = size(S,1);
    
    sbid = S(:,Vex(irxn))<0;    
    prid = S(:,Vex(irxn))>0;
    
    % remove protons
    sbid([he hc]) = 0;
    prid([he hc]) = 0;
    
    Kscol = zeros(nmet,1);
    Kpcol = zeros(nmet,1); 
    
    pvec = estimateKm(pvec,sbid,prid,mc,Kscol,Kpcol,Vex(irxn));
end

% set irreversible kcats
pvec.krev(~model.rev) = 0;
% for irxn = 1:length(model.rxns)
%     if ~model.rev(irxn)
%         pvec.krev(irxn)=0;
%     end
% end

% exhcnage reactions
pvec.kfwd(model.VFex) = 0;
pvec.krev(model.VFex) = 0;

% restore Vmax from backup
pvec.Vmax = newp.Vmax;

% estimate Vmax
if all(check(Vind)>0)
    pvec.Vmax(pvec.delGr==0) = 0;
    
    % for Vind
    % simple vmax = vss/ck
    for irxn = 1:length(Vind)
        if ~isnan(newp.Vmax(Vind(irxn)))
            pvec.Vmax(Vind(irxn)) = newp.Vmax(Vind(irxn));
        else
            [~,ck] = CKinetics(model,pvec,mc,Vind(irxn));
            if ck
                pvec.Vmax(Vind(irxn)) = model.Vss(Vind(irxn))/(3600*ck);
            else
                pvec.Vmax(Vind(irxn)) = 1;
            end
        end
    end
    
    % other reactions
    for irxn = 1:length(Vex)
        if ~isnan(newp.Vmax(Vex(irxn)))
            pvec.Vmax(Vex(irxn)) = newp.Vmax(Vex(irxn));
        else
            [~,tk] = TKinetics(model,pvec,mc,Vex(irxn));
            if tk
                pvec.Vmax(Vex(irxn)) = model.Vss(Vex(irxn))/(3600*tk);
            else
                pvec.Vmax(Vex(irxn)) = 1;
            end
        end
    end 
    % pvec.Vmax(Vex(~isnan(newp.Vmax(Vex)))) =...
    % newp.Vmax(Vex(~isnan(newp.Vmax(Vex))));
    
    % calculate fluxes for ETC reactions
    [~,etck] = ETCflux(model,pvec,mc,flux);
    rnlst = {'ATPS4r','NADH16','CYTBD'};
    for irxn = 1:length(rnlst)
        tfr = strcmpi(model.rxns,rnlst{irxn});
        if ~isnan(newp.Vmax(tfr))
            pvec.Vmax(tfr) = newp.Vmax(tfr);
        else
            if any(tfr) && etck(tfr)
                pvec.Vmax(tfr) = model.Vss(tfr)/(3600*etck(tfr));
            else
                pvec.Vmax(tfr) = 1;
            end
        end
    end    
       
    % atp maintanance
    tatpm = strcmpi(model.rxns,'ATPM');
    if any(tatpm)
        sbid = model.S(:,tatpm)<0;
        sbid(h2o) = 0;
        atpmk = 18.84*mc(sbid)/pvec.K(sbid,tatpm)/(1+mc(sbid)/pvec.K(sbid,tatpm));
        if logical(atpmk)
            pvec.Vmax(tatpm) = model.Vss(tatpm)/(3600*atpmk);
        else
            pvec.Vmax(tatpm) = 1;
        end
    end
    
    pvec.Vmax(model.VFex) = 0;    
    pvec.Vmax(model.Vss==0) = 0;
    pvec.feasible = 1;   
else
    fprintf('Thermodynamically infeasible parameters\n');
    fprintf('Discontinuing\n');
    pvec.feasible = 0;
    return
end

% check - calculate initial flux
flux = iflux(model,pvec,mc);

function check = checkthermo(fhandle,check)

%#check for vss and delGr direction      
flux(Vind(irxn)) = fhandle(model,pvec,mc,Vind(irxn));
%     flux = ETCflux(model,mc,flux);
if pvec.delGr(Vind(irxn)) ~= 0
    if flux(Vind(irxn))*pvec.delGr(Vind(irxn))<0
        check(Vind(irxn)) = 1;
    else
        check(Vind(irxn)) = -1;
    end    
else
    if flux(Vind(irxn))*pvec.delGr(Vind(irxn))<1e-6
        check(Vind(irxn)) = 1;
    else
        check(Vind(irxn)) = -1;
    end
end 

end

end

    
    
            



        