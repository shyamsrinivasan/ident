function newK = buildmodels(model,pvec,mc,rxn_add,rxn_excep,nmodels)
if nargin<6
    nmodels = 1;
end
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

% ETC rxns list
etcrxns = {'ATPS4r','NADH16','CYTBD','SUCDi','FRD7'};
etcrxns = cellfun(@(x)strcmpi(x,model.rxns),etcrxns,'UniformOutput',false);
etcrxns = cellfun(@(x)find(x),etcrxns,'UniformOutput',false);
etcrxns = cell2mat(etcrxns);

% parameters from backup in pvec(input)
delGr = pvec.delGr;
KIact = pvec.KIact;
KIihb = pvec.KIihb;
Vmax = pvec.Vmax;

nrxn = model.nt_rxn;
flux = zeros(nrxn,nmodels);

% sample nmodel Kms for all reactions in Vind
newK = samplesigma(model,mc,pvec.K,pvec.kfwd,pvec.krev,Vind,nmodels);

% other reactions 
% transport reactions x[e] <==> x[c]
new_excep = union(rxn_excep,model.rxns(Vind));
Vex = addToVind(model.rxns,model.Vex,[],new_excep);

% sample nmodel Kms for all reactions in Vex
newKVex = samplesigma(model,mc,pvec.K,pvec.kfwd,pvec.krev,Vex,nmodels);

% assign parameters to nmodels structure array parallel
if nmodels>10
    parfor im = 1:nmodels
        % set Vex rxn parameters in newK from newKVex
        newK(im).K(:,Vex) = newKVex(im).K(:,Vex);
        newK(im).kfwd(Vex) = newKVex(im).kfwd(Vex);
        newK(im).krev(Vex) = newKVex(im).krev(Vex);
        
        newK(im).delGr = delGr;
        newK(im).KIact = KIact;
        newK(im).KIihb = KIihb;  
        % check for consistency with fluxes in vss and delGr
        newK(im).check = checkthermo(model,newK(im),mc,Vind,@CKinetics); 
        
        % echange rxn flxs
        newK(im).kfwd(model.VFex) = 0;
        newK(im).krev(model.VFex) = 0;
        
        % determine Vmax values for all Vind rxns
        [newVmax,feasible] = getVmax(Vmax,model,newK(im),mc,Vind,@CKinetics);
        newK(im).Vmax(Vind) = newVmax(Vind);           
        
        % determine Vmax values for all Vex rxns
        newVmax = getVmax(Vmax,model,newK(im),mc,Vex,@TKinetics);
        newK(im).Vmax(Vex) = newVmax(Vex);    
        
        % determine Vmax values for other reactions in ETC
%         [newVmax,feasible] = getVmax(Vmax,model,newK(im),mc,Vex,@ETCflux);
        newK(im).feasible = feasible;  
    end
else
    for im = 1:nmodels
        % set Vex rxn parameters in newK from newKVex
        newK(im).K(:,Vex) = newKVex(im).K(:,Vex);
        newK(im).kfwd(Vex) = newKVex(im).kfwd(Vex);
        newK(im).krev(Vex) = newKVex(im).krev(Vex);
        
        newK(im).delGr = delGr;
        newK(im).KIact = KIact;
        newK(im).KIihb = KIihb;  
        % check for consistency with fluxes in vss and delGr
        newK(im).check = checkthermo(model,newK(im),mc,Vind,@CKinetics); 
        
        % echange rxn flxs
        newK(im).kfwd(model.VFex) = 0;
        newK(im).krev(model.VFex) = 0;
        
        % determine Vmax values for all Vind rxns
        [newVmax,feasible] = getVmax(Vmax,model,newK(im),mc,Vind,@CKinetics);
        newK(im).Vmax(Vind) = newVmax(Vind);           
        
        % determine Vmax values for all Vex rxns
        newVmax = getVmax(Vmax,model,newK(im),mc,Vex,@TKinetics);
        newK(im).Vmax(Vex) = newVmax(Vex);    
        
        % determine Vmax values for other reactions in ETC
%         [newVmax,feasible] = getVmax(Vmax,model,newK(im),mc,Vex,@ETCflux);
        newK(im).feasible = feasible;    
    
        % check - calculate initial flux
        flux(:,im) = iflux(model,newK(im),mc);
    end
end


        
        
 % calculate fluxes for ETC reactions
%     [~,etck] = ETCflux(model,pvec,mc,flux);
%     rnlst = {'ATPS4r','NADH16','CYTBD'};
%     for irxn = 1:length(rnlst)
%         tfr = strcmpi(model.rxns,rnlst{irxn});
%         if ~isnan(newp.Vmax(tfr))
%             pvec.Vmax(tfr) = newp.Vmax(tfr);
%         else
%             if any(tfr) && etck(tfr)
%                 pvec.Vmax(tfr) = model.Vss(tfr)/(3600*etck(tfr));
%             else
%                 pvec.Vmax(tfr) = 1;
%             end
%         end
%     end    


% sample nmodel kcat using Brigg's Haldane for all reactions in Vind
% pvec = samplekcat(model,pvec,sbid,prid,Vind(irxn),mc,kfwdbkup,kbkwbkup,rerun);

% for irxn = 1:nrxn
%     sbid = S(:,Vind(irxn))<0;    
%     prid = S(:,Vind(irxn))>0;        
%     % remove water
%     sbid(h2o) = 0;
%     prid(h2o) = 0;      
%     % remove protons
%     sbid([he hc]) = 0;
%     prid([he hc]) = 0;
%     
%     % no parameters for cofactors - assumed abundant 
%     % cofactros are assumed as compensated species
%     % hence   
%     rerun = 0;
%     Ksbackup = pvec.K(sbid,Vind(irxn));
%     Kpbackup = pvec.K(prid,Vind(irxn));
%     kfwdbkup = pvec.kfwd(Vind(irxn));
%     kbkwbkup = pvec.krev(Vind(irxn));
%     while check(Vind(irxn))~=1 
% %         if any(Ksbackup==1)||any(Kpbackup==1)
%             % sampling of parameters needs to be recursive until check (see below) is 1
%         pvec = estimateKm(pvec,sbid,prid,mc,Ksbackup,Kpbackup,Vind(irxn),rerun);
%     
%     % forward and backward catalytic rates
%     % kfwd and kbkw
%     % kfwd or kbkw is sampled basedon reaction directionality from FBA for
%     % thermodynamic consistency
%     % sampling done only for unknown values
% %     fprintf('%s \t delG = %3.6g \t Vflux = %3.6g\t',model.rxns{Vind(irxn)},...
% %              pvec.delGr(Vind(irxn)),...
% %              model.Vss(Vind(irxn)));  
% %     if ~any(strcmpi(model.rxns{Vind(irxn)},'ATPS4r'))
%         % ATPsynthase does not strictly obey Briggs Haldane
%         pvec = samplekcat(model,pvec,sbid,prid,Vind(irxn),mc,kfwdbkup,kbkwbkup,rerun);
% %     end    
%         pvec.Vmax(Vind(irxn)) = 1;
%         pvec.Vmax(model.Vss==0) = 0;
%     
%         % check for vss and delGr direction
%         check = checkthermo(@CKinetics,check);
%         rerun = 1;
%         % if not satisfied => resample -> use a while loop?
% %         else
% %             fprintf('Parameter estimation entering an infinite loop for %s\n',model.rxns(Vind(irxn)))
% %         end
%     end
% end



% for irxn = 1:length(Vex)

%     
%     pvec = estimateKm(pvec,sbid,prid,mc,Kscol,Kpcol,Vex(irxn));
% end






% estimate Vmax
% if all(check(Vind)>0)
%        
%        
%     % atp maintanance
%     tatpm = strcmpi(model.rxns,'ATPM');
%     if any(tatpm)
%         sbid = model.S(:,tatpm)<0;
%         sbid(h2o) = 0;
%         atpmk = 18.84*mc(sbid)/pvec.K(sbid,tatpm)/(1+mc(sbid)/pvec.K(sbid,tatpm));
%         if logical(atpmk)
%             pvec.Vmax(tatpm) = model.Vss(tatpm)/(3600*atpmk);
%         else
%             pvec.Vmax(tatpm) = 1;
%         end
%     end
%     
%     pvec.Vmax(model.VFex) = 0;    
%     pvec.Vmax(model.Vss==0) = 0;
%     pvec.feasible = 1;   
% else
%     fprintf('Thermodynamically infeasible parameters\n');
%     fprintf('Discontinuing\n');
%     pvec.feasible = 0;
%     return
% end





end

    
    
            



        