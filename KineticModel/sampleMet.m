function [MC_sampl] = sampleMet(FBAmodel,ism)
nmetab = FBAmodel.nt_metab;
pd = makedist('Beta');
MClow = FBAmodel.MClow(1:nmetab);
MChigh = FBAmodel.MChigh(1:nmetab);
mSample = MClow + (MChigh - MClow).*random(pd,nmetab,1);
MC_sampl = zeros(nmetab,1);
MC_sampl = mSample;

% MC_sampl(MClow ~= MChigh) = mSample(MClow ~= MChigh);
% %ATP AMP ADP
% ec = 0.8;
% ATP = strcmpi('atp[c]',FBAmodel.mets);
% ADP = strcmpi('adp[c]',FBAmodel.mets);
% AMP = strcmpi('amp[c]',FBAmodel.mets);
% MC_sampl(ATP) = (MC_sampl(AMP)-MC_sampl(ADP)*(1/(2*ec)-1))/(1/ec-1);

%Sample mets based on Lsawrence's Noisy Metabolomics Sampling
%Methodology (Directly load a pre-sampled file)
% ism = 10;
% %Model
% load('C:\Users\shyam\Documents\Courses\CHE1125Project\mat_files\KineticModel\Samples\centralIshiiModel.mat');
% %Samples
% load('C:\Users\shyam\Documents\Courses\CHE1125Project\mat_files\KineticModel\Samples\CentralIshii_1.mat');
% mets = cellfun(@changeN,lower(model.mets),'UniformOutput',0);
% MCsample = zeros(length(FBAmodel.mets),1000);
% for im = 1:length(FBAmodel.mets)
%     tfm = strcmpi(FBAmodel.mets{im},mets);
%     if any(tfm)
%         fprintf('%d Metabolite:%s\n',im,FBAmodel.mets{im});
%         MCsample(im,:) = exp(points(tfm,:));%points from loading mat file #2
%         MC_sampl(im) = MCsample(im,ism);  
%     else
%         MCsample(im,:) = MClow(im) + (MChigh(im) - MClow(im)).*random(pd,1,1000);
%         MC_sampl(im) = MCsample(im,ism); 
%     end    
% end
H2o = find(strcmpi('h2o[c]',FBAmodel.mets));
Hion = find(strcmpi('h[c]',FBAmodel.mets));
H2oe = find(strcmpi('h2o[e]',FBAmodel.mets));
Hione = find(strcmpi('h[e]',FBAmodel.mets));

nad = [];%strcmpi('nad[c]',FBAmodel.mets);
nadh = [];%strcmpi('nadh[c]',FBAmodel.mets);

if any(MClow == MChigh)
    MC_sampl(MClow == MChigh) = MClow(MClow == MChigh);   
end
if ~isempty(H2oe)
    MC_sampl(H2oe) = 1000;
end

MC_sampl(nad) = 1000*MC_sampl(nadh);
% MC = MC_sampl;
% [nm,nr] = size(FBAmodel.S);
% for irxn = 1:nr
%     if irxn ~= [12,48]
%     Keq = FBAmodel.Keq(irxn);
%     Vss = FBAmodel.Vss(irxn);
%     subs = FBAmodel.S(:,irxn)<0;
%     prud = FBAmodel.S(:,irxn)>0;
%     subs([H2o Hion H2oe Hione]) = 0;
%     prud([H2o Hion H2oe Hione]) = 0;
%     if any(subs) && any(prud)
%         MCsub = MC(subs);       
%         MCprd = MC(prud);        
%         gamma = prod(MCsub)-prod(MCprd)/Keq;
%         if Vss ~= 0
%             while gamma*Vss < 0 
%                 MCsub = FBAmodel.MClow(subs) +...
%                         (FBAmodel.MChigh(subs) - FBAmodel.MClow(subs)).*...
%                         betarnd(1.5,4.5,length(find(subs)),1);
%                 MCprd = FBAmodel.MClow(prud) +...
%                         (FBAmodel.MChigh(prud) - FBAmodel.MClow(prud)).*...
%                         betarnd(1.5,4.5,length(find(prud)),1);
%                 gamma = prod(MCsub)-prod(MCprd)/Keq;
%             end
%         end
%         MC(subs) = MCsub;
%         MC(prud) = MCprd;
%     end    
%     end    
% end
% for irxn = 1:length(FBAmodel.Vind)
%     subs = FBAmodel.S(:,FBAmodel.Vind(irxn))<0;
%     prud = FBAmodel.S(:,FBAmodel.Vind(irxn))>0;
%     if any(subs) && any(prud)
%         MCsub = MC(subs);
%         Msm_sub = MCsample(subs,ism);
%         MCprd = MC(prud);
%         Msm_prd = MCsample(prud,ism);
%         if any(MCsub==0) && ~any(MCprd==0)
%             MCsub(MCsub==0) = prod(MCprd)/FBAmodel.Keq(FBAmodel.Vind(irxn));
%             MC(subs) = MCsub;
%         end  
%         if ~any(MCsub==0) && any(MCprd==0)
%             MCprd(MCprd==0) = prod(MCsub)*FBAmodel.Keq(FBAmodel.Vind(irxn));
%             MC(prud) = MCprd;
%         end
%         if any(MCsub==0) && any(MCprd==0)
%             MCsub(MCsub==0) = Msm_sub(MCsub==0);
%             MCprd(MCprd==0) = prod(MCsub)*FBAmodel.Keq(FBAmodel.Vind(irxn));
%             MC(subs) = MCsub;
%             MC(prud) = MCprd;
%         end        
%     end
% end
% MC_sampl = MC;
% [nm,nr] = size(model.S);
% [nm2,nr2] = size(FBAmodel.S);
% for ir = 1:nr
%     metid = mets(logical(model.S(:,ir)));
%     display(metid);
% end
% newRxn = cell(nr2,1);
% oldRxn = cell(nr2,1);
% for ir = 1:nr
%     metid = mets(logical(model.S(:,ir)));
%     for irxn = 1:nr2
%         metid2 = FBAmodel.mets(logical(FBAmodel.S(:,irxn)));  
%         if length(metid2) == length(metid)
%         met_diff = setdiff(metid2,metid);
%         if isempty(met_diff)
%             newRxn{irxn} = model.rxns{ir};
%             oldRxn{irxn} = FBAmodel.Enzyme{irxn};
%         end
%         end
%     end
% end


%ATP AMP ADP
ec = 0.8;
ATP = strcmpi('atp[c]',FBAmodel.mets);
ADP = strcmpi('adp[c]',FBAmodel.mets);
AMP = strcmpi('amp[c]',FBAmodel.mets);
% MCsample(ATP,:) = (MCsample(AMP,:)-MCsample(ADP,:).*(1/(2*ec)-1))./(1/ec-1);
% MC_sampl(ATP) = MCsample(ATP,ism);


% MC_sampl = [MC_sampl;zeros(model.next_metab,1)];
%Check for thermodynamic consistency
% for irxn = 1:model.n_rxn
%     
% end
%non Linear constraint optimzation to sample concentrations
% X = fmincon(0,zeros(nmetab,1),[],[],[],[],model.MClow,model.MChigh,@nonlncon);
    
return
function newN = changeN(str)
if ~isempty(regexp(str,'\w(?:xt|EXT|ext)$','ONCE'))
    newN = [str(1:regexp(str,'\w(?:xt|EXT|ext)$')) '[e]'];
else
    newN = [str '[c]'];
end
return
