function [MC_sampl,KVl] = sampleMet(FBAmodel,MC_sampl)

nint_metab = FBAmodel.nint_metab;
next_metab = FBAmodel.next_metab;
nmetab = FBAmodel.nt_metab;
nrxn = FBAmodel.nt_rxn;
pd = makedist('Beta','a',2,'b',2);
MClow = FBAmodel.MClow(1:nmetab);
MChigh = FBAmodel.MChigh(1:nmetab);

%sample only zero metabolites
zro_m = logical(~MC_sampl);
nzr_m = logical(MC_sampl);

nzro = length(find(zro_m));
%sample metabolites
mSample = MClow + (MChigh - MClow).*random(pd,nmetab,1);
% MC_sampl = mSample;
MC_sampl(zro_m) = mSample(zro_m);

%sample KVl - Velocity constant for CKM
low = 0.1;
high = 1000;
KVl = low + (high-low).*rand(nrxn,1);

% MC_sampl(MClow ~= MChigh) = mSample(MClow ~= MChigh);


%Sample mets based on Lsawrence's Noisy Metabolomics Sampling
%Methodology (Directly load a pre-sampled file)

ism = 100;
%Model
load('C:\Users\shyam\Documents\Courses\CHE1125Project\mat_files\KineticModel\Samples\centralIshiiModel.mat');
%Samples
load('C:\Users\shyam\Documents\Courses\CHE1125Project\mat_files\KineticModel\Samples\CentralIshii_1.mat');
mets = cellfun(@changeN,lower(model.mets),'UniformOutput',0);
MCsample = zeros(length(FBAmodel.mets),1000);
% for im = 1:length(FBAmodel.mets)
%     tfm = strcmpi(FBAmodel.mets{im},mets);
%     if any(tfm)
% %         fprintf('%d Metabolite:%s\n',im,FBAmodel.mets{im});
%         MCsample(im,:) = exp(points(tfm,:));%points from loading mat file #2
%         MC_sampl(im) = MCsample(im,ism);  
%     else
MCsample(zro_m,:) = repmat(MClow(zro_m),[1,1000]) +...
                   (repmat(MChigh(zro_m),[1,1000]) -...
                    repmat(MClow(zro_m),[1,1000])).*random(pd,nzro,1000);
MC_sampl(zro_m) = MCsample(zro_m,ism); 
%     end    
% end
H2o = find(strcmpi('h2o[c]',FBAmodel.mets));
Hion = find(strcmpi('h[c]',FBAmodel.mets));
H2oe = find(strcmpi('h2o[e]',FBAmodel.mets));
Hione = find(strcmpi('h[e]',FBAmodel.mets));

nad = [];%strcmpi('nad[c]',FBAmodel.mets);
nadh = [];%strcmpi('nadh[c]',FBAmodel.mets);

% if any(MClow == MChigh)
%     MC_sampl(MClow == MChigh) = MClow(MClow == MChigh);   
% end
if ~isempty(H2oe) && MC_sampl(H2oe) == 0
    MC_sampl(H2oe) = 1000;
end

MC_sampl(nad) = 1000*MC_sampl(nadh);

% %ATP AMP ADP
ec = 0.8;
ATP = strcmpi('atp[c]',FBAmodel.mets);
ADP = strcmpi('adp[c]',FBAmodel.mets);
AMP = strcmpi('amp[c]',FBAmodel.mets);
MC_sampl(ATP) = (MC_sampl(AMP)-MC_sampl(ADP)*(1/(2*ec)-1))/(1/ec-1);


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
