function [MC_sampl,KVl] = sampleMet(FBAmodel,ism)

nmetab = FBAmodel.nt_metab;
nrxn = FBAmodel.nt_rxn;
pd = makedist('Beta','a',2,'b',2);
MClow = FBAmodel.MClow(1:nmetab);
MChigh = FBAmodel.MChigh(1:nmetab);

%sample metabolites
mSample = MClow + (MChigh - MClow).*random(pd,nmetab,1);
% MC_sampl = zeros(nmetab,1);
MC_sampl = mSample;

%sample KVl - Velocity constant for CKM
low = 0.1;
high = 1000;
KVl = low + (high-low).*rand(nrxn,1);

% MC_sampl(MClow ~= MChigh) = mSample(MClow ~= MChigh);


%Sample mets based on Lsawrence's Noisy Metabolomics Sampling
%Methodology (Directly load a pre-sampled file)

ism = 10;
%Model
load('C:\Users\shyam\Documents\Courses\CHE1125Project\mat_files\KineticModel\Samples\centralIshiiModel.mat');
%Samples
load('C:\Users\shyam\Documents\Courses\CHE1125Project\mat_files\KineticModel\Samples\CentralIshii_1.mat');
mets = cellfun(@changeN,lower(model.mets),'UniformOutput',0);
MCsample = zeros(length(FBAmodel.mets),1000);
for im = 1:length(FBAmodel.mets)
%     tfm = strcmpi(FBAmodel.mets{im},mets);
%     if any(tfm)
% %         fprintf('%d Metabolite:%s\n',im,FBAmodel.mets{im});
%         MCsample(im,:) = exp(points(tfm,:));%points from loading mat file #2
%         MC_sampl(im) = MCsample(im,ism);  
%     else
        MCsample(im,:) = MClow(im) + (MChigh(im) - MClow(im)).*random(pd,1,1000);
        MC_sampl(im) = MCsample(im,ism); 
%     end    
end
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
