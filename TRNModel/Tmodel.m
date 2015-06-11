% function [trnmodel,defparval,ngene,ngap,tnreg,regprotein,rgene,rmetab,rprot,...
%           defvalueGene,C] = Tmodel(fname1,fname2,defparval)
%**************************************************************************
%Model construction from text file for trnascriptional regulatory networks.
%All coefficients and rates are assigned default values if not provided in
%file
%*******Input(s)
%fname1 - 
%fname2 - 
%*******Optional Input(s)
%defparval - Structure of default parameter values to be used
%*******Output(s)
%trnmodel - Matlab structure consisting of the transcriptional regulatory
%model from text file
%defparval - 
%ngene - Total number of genes in model
%ngap - Total number of non-gene associate proteins (Proteins associated
%with metaolites - simplification of the TCS for regulation in bacteria
%tnreg - Total number of regulators
%Rest to be completed
%December 12th 2013
%Extracting Activation & Repression Coefficients in One go
%Need for if nargin < x arguments within the function

%September 2013 - version 1.0
%December 2013 - Changed how srate is represented on a per gene basis as
%opposed to a matrix form corresponding to trnmodel.RS
%December 19th 2013 - Reverted srate back to TF basis
%August 2014 - version 2.1
%**************************************************************************
function [model,FBAmodel,defparval,C] =...
         Tmodel(fname1,FBAmodel,pmeter,variable,defparval)

if nargin < 5
    %Default parameters
    defparval.accoeff = 10;%mmole
    defparval.repcoeff = 10;%mmole
    defparval.rephill = 2;
    defparval.dualcoeff = 1e-3;%mmole
    defparval.brate = 1.66e+9;
    defparval.trate = 1;%trate is growth dependent => connection matrix
    defparval.ptrate = 0.06;
    defparval.mdecay = 0.0031;%s-1
    defparval.pdecay = 3.83e-6;%s-1    
else
    if ~isfield(trate)
        defparval.trate = 1;%s-1
    end
    if ~isfield(mdecay)
       defparval.mdecay = 0.0031;%s-1
    end
    if ~isfield(pdecay)
        defparval.pdecay = 3.83e-6;%s-1    
    end    
end
     
if nargin < 1 
    fprintf('Atleast 1 file name Missing.\n Cannot Proceed\n');    
    return
end    
fileid1 = fopen(fname1);
% fileid2 = fopen(fname2);
if fileid1 == -1
    fprintf('File %s cannot be opened.', fname1);
    model = struct([]);
    return;
end
C = textscan(fileid1, '%s%s%s%s%s%s%f%f%f%f%f%s%s%s%s', 'Delimiter', '\t',...
            'TreatAsEmpty', {'None'}, 'HeaderLines', 1);
fclose(fileid1);
% fclose(fileid2);
model.Gene = C{1};
%ecocyc = C{2};
prot_name = C{3};
gene_prot = C{4};
gene_metab = C{5};
p_inducer = C{12};
m_inducer = C{13};

nt_gene = length(model.Gene);
nt_metab = FBAmodel.nt_metab;
nint_metab = FBAmodel.nint_metab;
next_metab = FBAmodel.next_metab;
nt_prot = length(FBAmodel.Enzyme)+1;%+1 for RNAP
model.Enzyme = [FBAmodel.Enzyme;'RNAp'];
model.Metabolites = FBAmodel.Metabolites;
bm_ind = strcmpi(model.Metabolites,'Biomass');
%Regulators = [Proteins;Metabolites w/o biomass];
if any(bm_ind) && (length(model.Metabolites) > length(bm_ind))
    model.Regulators = unique([model.Enzyme;...
                               FBAmodel.Metabolites(1:find(bm_ind)-1);...
                               FBAmodel.Metabolites(find(bm_ind)+1:end)]);
elseif any(bm_ind)
    model.Regulators = unique([model.Enzyme;...
                               FBAmodel.Metabolites(1:find(bm_ind)-1)]);
else
    model.Regulators = unique([model.Enzyme;FBAmodel.Metabolites]);
end
nt_reg = length(model.Regulators);
model.GeneRules = cell(nt_gene,1);
model.RS = sparse(nt_gene,nt_reg);
model.Coefficient = sparse(nt_gene,nt_reg);
model.trate = sparse(nt_reg,nt_gene);
model.brate = zeros(nt_gene,1);
model.Kb = zeros(nt_reg,1);
model.Kub = zeros(nt_reg,1);
%Steady state [mRNA]
model.SSmRNA = zeros(nt_gene,1);
%Metabolite Concentraton
if isfield(variable,'MC')
    model.MC = zeros(nt_metab,1);
    model.MC(1:nint_metab) = variable.MC;
%     model.MC(nint_metab+1:nint_metab+next_metab) = FBAmodel.M*variable.MC;
else
    model.MClow = FBAmodel.MClow;
    model.MChigh = FBAmodel.MChigh;
end
%Protein Concentration
model.SSprot = zeros(nt_prot,1);
model.pmRatio = zeros(nt_prot,1);
for iprot = 1:nt_prot   
    tf_prot = strcmpi(model.Enzyme{iprot},prot_name);
    if any(tf_prot)
        %Initial Steady state protein concentration for normalization
        if ~isempty(C{11}(tf_prot)) && ~isnan(C{11}(tf_prot))
            model.SSprot(iprot) = C{11}(tf_prot);
        end
        %Protein/mRNA ratio to calculate [mRNA] at Steady state
        if ~isempty(C{10}(tf_prot)) && ~isnan(C{10}(tf_prot))
            model.pmRatio(iprot) = C{10}(tf_prot);
        end
    end
end
%RNAP
model.SSprot(end) = model.SSprot(iprot-1);%umole
model.pmRatio(end) = 1000;

coeff_flag = 0;
for igene = 1:nt_gene
    assign_coeff = 0;
    %Assign Protein to Gene    
    ireg = length(model.Regulators)+1;
    iprot = length(model.Enzyme)+1;
    tfg_all = strcmpi(prot_name{igene},model.Regulators);    
    if ~any(tfg_all) 
        model.Regulators{ireg} = prot_name{igene};
        model.RS(igene,ireg) = 0;
        model.Coefficient(igene,ireg) = 0;
        model.trate(ireg,igene) = defparval.trate; 
    else
        model.trate(tfg_all,igene) = defparval.trate;
    end
    tfgp = strcmpi(prot_name{igene},model.Enzyme);
    if ~any(tfgp)
        model.Enzyme{iprot} = prot_name{igene};
    end       
    %Coefficients    
    if ~isempty(C{6}{igene})
        [coeff] = regexp(C{6}{igene},'(\w*\w+.?\S\d+)\((\W+.?)\)+','tokens');
        assign_coeff = 1;
    else
        coeff = {};
    end    
    %Binding/Unbinding Coefficients for Complexes
    if ~isempty(p_inducer{igene}) || ~isempty(m_inducer{igene})
        if ~isempty(C{14}{igene})
            kbind = strsplit(C{14}{igene},',');                        
        else
            kbind = {};
        end
        if ~isempty(C{15}{igene})
            kubind = strsplit(C{15}{igene},',');            
        else
            kubind = {};
        end
    end
    %Build Regulatory Rules           
    ruleterm = {}; 
    tterms = 0;    
    %Proteins
    if ~isempty(gene_prot{igene}) &&...
        (~isempty(p_inducer{igene})||~isempty(m_inducer{igene}))  
        %Protein Complexes
        [ruleterm,tterms] =...
        form_rule(ruleterm,gene_prot{igene},length(ruleterm),1,tterms,...
                  p_inducer{igene},m_inducer{igene},kbind,kubind);
    else                          
        [ruleterm,tterms] =...
        form_rule(ruleterm,gene_prot{igene},length(ruleterm),1,tterms);
    end
    %Metabolites 
    if ~isempty(gene_metab{igene}) &&...
        (~isempty(p_inducer{igene})||~isempty(m_inducer{igene}))       
        %Metabolite Complexes
        [ruleterm,tterms] =...
        form_rule(ruleterm,gene_metab{igene},length(ruleterm),2,tterms,...
                  p_inducer{igene},m_inducer{igene},kbind,kubind);
    else
        [ruleterm,tterms] =...
        form_rule(ruleterm,gene_metab{igene},length(ruleterm),2,tterms);
    end
    model.GeneRules{igene} = ruleterm;   
    %Basal Transcription Rate - b in (1/b+pact);
    if ~isempty(C{7}(igene))
        if ~isnan(C{7}(igene))
            model.brate(igene) = C{7}(igene);        
        else
            model.brate(igene) = defparval.brate;      
        end
    end    
    %mRNA Concentrations
    if ~isempty(C{8}(igene))
        if ~isnan(C{8}(igene))
            %C{8}(igene);
            model.SSmRNA(igene) = model.SSprot(igene)/model.pmRatio(igene);            
        end
    end   
end

model.S = FBAmodel.S;
model.SI = FBAmodel.SI;
model.MClow = FBAmodel.MClow;
model.MChigh = FBAmodel.MChigh;
model.Vss = FBAmodel.Vss;
model.Kcat = FBAmodel.Kcat;
model.delSGr = FBAmodel.delSGr;
model.Keq = FBAmodel.Keq;
model.pmeter = pmeter;
%Adding active/Inactive states for complexes (non regulating - Eg. LacI and
%LacI* in LacI-lac[c]
noregcplx = model.Metabolites(~cellfun('isempty',strfind(model.Metabolites,'-')));
n_noregclx = length(noregcplx);
for icpx = 1:n_noregclx
    %Add corresponding active state
    regstr = strsplit(noregcplx{icpx},'-');    
    %regstr{1} - Regulator
    tfr_c = strcmpi(regstr{1},model.Regulators);
    %regstr{2} - Deregulating Operator
    if length(regstr) > 1
    tfr_o = strcmpi(regstr{2},model.Regulators);    
        if any(tfr_c) && any(tfr_o)
            %Active form of regulator
            act_reg = sprintf('%s*',regstr{1});
            %Replace regstr{1} with act_reg and make changes to model
            model = replaceRegulator(model,regstr{1},act_reg);
        end 
    end
end

nt_reg = length(model.Regulators);
%model.Reg Assignment
%remove common elements from metrabolites and proteins to identify indices
%prot_ind and metab_ind correctly below
metab_ind = zeros(nt_reg,1);
prot_ind = zeros(nt_reg,1);
cm_ind = zeros(nt_reg,1);
cplx_ind = zeros(nt_reg,1);
model.SSreg = zeros(nt_reg,1);

for kreg = 1:nt_reg
    prot = 0;metab = 0;
    if any(strcmpi(model.Regulators{kreg},model.Metabolites))
        metab = 1;
    end
    if any(strcmpi(model.Regulators{kreg},model.Enzyme))
       prot = 1;
    end
    if prot && metab%Common to both 
%         if isfield(variable,'MC')
%             if any(bm_ind)
%                 tfmc = strcmpi(model.Regulators{kreg},...
%                                [model.Metabolites(1:find(bm_ind)-1),...
%                                 model.Metabolites(find(bm_ind)+1:end)]);
%             else
%                 tfmc = strcmpi(model.Regulators{kreg},model.Metabolites);
%             end
%             model.SSreg(kreg) = model.MC(tfmc);            
%         else
            model.SSreg(kreg) =...
            model.SSprot(strcmpi(model.Regulators{kreg},model.Enzyme));
%         end
        cm_ind(kreg) = 1;
    elseif prot && ~metab
        model.SSreg(kreg) =...
        model.SSprot(strcmpi(model.Regulators{kreg},model.Enzyme));
        prot_ind(kreg) = 1;
    elseif ~prot && metab
        %Non Regulating Complexes - LacI-Allolactose
%         noRegcmp = strsplit(model.Regulators{kreg},'-');
%         if length(noRegcmp)>1%If Regulator{kreg} is a complex
%             model.SSreg(kreg) = 0;
%         else
            %Regular Metabolites
            if isfield(variable,'MC')
                 if ~isempty(bm_ind)
                    tfmc = strcmpi(model.Regulators{kreg},...
                                   [model.Metabolites(1:bm_ind-1),...
                                    model.Metabolites(bm_ind+1:end)]);
                else
                    tfmc = strcmpi(model.Regulators{kreg},model.Metabolites);
                end
                model.SSreg(kreg) = model.MC(tfmc);    
            end
%         end
        metab_ind(kreg) = 1;
    elseif ~prot && ~metab%Not in Both Regulatory Complexes
        model.SSreg(kreg) = 0;  
        cplx_ind(kreg) = 1;
    end
end

prot_ind = logical(prot_ind);
metab_ind = logical(metab_ind);
cm_ind = logical(cm_ind);
cplx_ind = logical(cplx_ind);


% model.M = FBAmodel.M;

%Change Enzyme Order
[~,bmrxn] = find(model.S(strcmpi(model.Metabolites,'Biomass'),:) > 0);
if any(bmrxn)
    bm_prot = strcmpi(model.Enzyme{bmrxn},model.Regulators);
else
    bm_prot = [];
end
rnape_ind = strcmpi('RNAp',model.Regulators);
prot_ind(rnape_ind) = 0;
prot_ind(bm_prot) = 0;
newEnzyme = [model.Regulators(prot_ind);...
             model.Regulators(bm_prot);...
             model.Regulators(cm_ind);...
             model.Regulators(rnape_ind)];

newEindx = zeros(length(newEnzyme),1);
for inewE = 1:length(newEnzyme)
    tf_newE = strcmp(newEnzyme{inewE},model.Enzyme);
    if any(tf_newE)
        newEindx(inewE) = find(tf_newE);
    end
end
rnap_ind = find(strcmpi(newEnzyme{inewE},model.Enzyme));
newEindx(newEindx==rnap_ind)=0;
newEindx(newEindx==0)=[];
model.Enzyme = newEnzyme;
model.Vss = [model.Vss(newEindx);zeros(length(rnap_ind),1)];
model.Kcat = [model.Kcat(newEindx);zeros(length(rnap_ind),1)];
model.delSGr = [model.delSGr(newEindx);zeros(length(rnap_ind),1)];
model.Keq = [model.Keq(newEindx);zeros(length(rnap_ind),1)];

%Adjust columns of S, SI, K and KI
% nt_rxn = size(model.S,2);
nt_rxn = FBAmodel.nt_rxn;
nt_m = size(model.S,1);
newS = [model.S(:,newEindx) sparse(nt_m,1)];
newK = [model.pmeter.K(:,newEindx) sparse(nt_m,1)];
newSI = [model.SI(:,newEindx) sparse(nt_m,1)];
newKIact = [model.pmeter.KIact(:,newEindx) sparse(nt_m,1)];
newKIihb = [model.pmeter.KIihb(:,newEindx) sparse(nt_m,1)];
nt_enz = size(newS,2);
%Chnage Metabolite Order       
bm_ind = strcmpi(model.Metabolites,'Biomass');
%External Metabolite Indices
exind_M = ~cellfun('isempty',regexp(model.Regulators,'\w(?:\[e\])$'));
% exind_M = ~cellfun('isempty',regexp(model.Regulators,'\w(?:xt)$'));
metab_ind(exind_M) = 0;
newMetab = [model.Regulators(cm_ind);...
            model.Regulators(metab_ind);...
            model.Regulators(exind_M);...
            model.Metabolites(bm_ind)];  

%Rearrange Metabolites
newMindx = zeros(length(newMetab),1);        
for inewM = 1:length(newMetab)
    tf_newm = strcmpi(newMetab{inewM},model.Metabolites);
    if any(tf_newm)
        newMindx(inewM) = find(tf_newm);
    end
end
% metabolites = model.Metabolites;
model.Metabolites = newMetab;
newS = newS(newMindx,:);
newK = newK(newMindx,:);
pmeter.K = newK;
newSI = newSI(newMindx,:);
%Activating
newKIact = newKIact(newMindx,:);
pmeter.KIact = newKIact;
%Inhibiting
newKIihb = newKIihb(newMindx,:);
pmeter.KIihb = newKIihb;
%Concentrations
model.MC = model.MC(newMindx,1);
       
model.RS = [model.RS(:,prot_ind),...
            model.RS(:,bm_prot),...
            model.RS(:,cplx_ind),...
            model.RS(:,cm_ind),... 
            model.RS(:,rnape_ind),...
            model.RS(:,metab_ind),...            
            model.RS(:,exind_M)];
model.Coefficient = [model.Coefficient(:,prot_ind),...                     
                     model.Coefficient(:,bm_prot),...
                     model.Coefficient(:,cplx_ind),...
                     model.Coefficient(:,cm_ind),...
                     model.Coefficient(:,rnape_ind),...   
                     model.Coefficient(:,metab_ind),...                     
                     model.Coefficient(:,exind_M)];
model.Regulators = [model.Regulators(prot_ind);...
                    model.Regulators(bm_prot);...
                    model.Regulators(cplx_ind);...
                    model.Regulators(cm_ind);... 
                    model.Regulators(rnape_ind);...
                    model.Regulators(metab_ind);...                    
                    model.Regulators(exind_M)];               
model.trate = [model.trate(prot_ind,:);...
               model.trate(bm_prot,:);...
               model.trate(cplx_ind,:);...
               model.trate(cm_ind,:);... 
               model.trate(rnape_ind,:);...
               model.trate(metab_ind,:);...               
               model.trate(exind_M,:)];
model.Kb =    [model.Kb(prot_ind,1);...
               model.Kb(bm_prot,1);...
               model.Kb(cplx_ind,1);...
               model.Kb(cm_ind,1);...
               model.Kb(rnape_ind,1);...
               model.Kb(metab_ind,1);...               
               model.Kb(exind_M,1)];
model.Kub =    [model.Kub(prot_ind,1);...
                model.Kub(bm_prot,1);...
                model.Kub(cplx_ind,1);...
                model.Kub(cm_ind,1);...
                model.Kub(rnape_ind,1);...
                model.Kub(metab_ind,1);...                
                model.Kub(exind_M,1)];           
model.SSreg = [model.SSreg(prot_ind,1);...
               model.SSreg(bm_prot,1);...  
               model.SSreg(cplx_ind,1);...
               model.SSreg(cm_ind,1);... 
               model.SSreg(rnape_ind,1);...
               model.SSreg(metab_ind,1);...               
               model.SSreg(exind_M)];
model.ptrate = model.trate;
model.ptrate(model.ptrate == 1) = defparval.ptrate;
model.S = newS;
model.SI = newSI;
%Recalculate all indices wrt new order of metabolites - Jan 26 2015
model.nt_rxn = size(model.S,2);
nt_rxn = model.nt_rxn;
model.nt_metab = length(model.Metabolites);
model.n_metab = length(find(cellfun('isempty',regexp(model.Metabolites,'\w(?:\[e\])$'))));
exind = ~cellfun('isempty',regexp(model.Metabolites,'\w(?:\[e\])$'));
model.bm_ind = find(strcmpi(model.Metabolites,'Biomass'));
Vuptake = [];
Vexind = [];
[~,allrxns] = find(newS(:,1:nt_enz));
% [~,allrxns] = find(newS(:,1:nt_rxn));
all_unqrxns = unique(allrxns);
nmetab_allrxns = histc(allrxns,all_unqrxns);
ex_rxn = all_unqrxns(nmetab_allrxns == 1);
[~,rxn] = find(newS(exind,ex_rxn));
Vext = ex_rxn(rxn);%one sided rxns using only external metabolites
[~,rxn1] = find(newS(~exind,ex_rxn));
ex_rxn = ex_rxn(rxn1);

for jrxn = 1:length(ex_rxn)
    if ~any(newS(:,ex_rxn(jrxn))<0)
        Vuptake = [Vuptake;ex_rxn(jrxn)];%one sided uptake rxn
    end
    if ~any(newS(:,ex_rxn(jrxn))>0)
        Vexind = [Vexind;ex_rxn(jrxn)];%one sided excretion rxn
    end
end
[~,rxn_enz] = find(newS);
Vd = setdiff(unique(rxn_enz),[Vuptake',Vexind',Vext']);
% Vd = setdiff(1:nt_rxn,[Vuptake',Vexind',Vext']);
[~,rxn] = find(newS(exind,Vd));
Vex = Vd(rxn);%two sided exchange reactions
[~,bmrxn] = find(newS(strcmpi(model.Metabolites,'Biomass'),:) > 0);
%find columns of newS which have no reactions
norxn = setdiff(1:nt_enz,allrxns);
%non Enzymatic reactions
noEnzreact = find(~cellfun('isempty',regexp(model.Enzyme,'(\w+.?)(?:_noenz)$')));
model.VnoEnz = noEnzreact;
%intracellular rxns
model.Vind = setdiff(1:nt_enz,[Vuptake;Vexind;bmrxn;Vext;Vex;norxn';noEnzreact]);
model.Vupind = Vuptake;
model.Vexind = Vexind;
model.Vex = Vex;
model.bmrxn = bmrxn;
model.rnap_ind = find(strcmpi('RNAp',model.Regulators));
[~,allactrxn] = find(newSI(:,1:nt_rxn)>0);
model.Vact_ind = unique(allactrxn);
[~,allihbrxn] = find(newSI(:,1:nt_rxn)<0);
model.Vihb_ind = unique(allihbrxn);
model.nt_rxn = nt_rxn;
model.n_rxn = length(model.Vind);
model.nt_reg = length(model.Regulators);

nprot_only = length(find(prot_ind))+length(find(bm_prot));
%Protein Complexes, Parent Proteins/Metabolites/Inducer Index - March 2015
%All complexes/algebraic variables
algind = ~cellfun('isempty',strfind(model.Regulators,'-'))|...
          ~cellfun('isempty',strfind(model.Regulators,'*'));
clxind = find(~cellfun('isempty',strfind(model.Regulators,'-')));      
model.CSA = sparse(length(model.Regulators),length(clxind));
for ic = 1:length(clxind)
    prot_cx = strsplit(model.Regulators{clxind(ic)},'-');
    tf_A = strcmpi(prot_cx{1},model.Regulators);
    if any(tf_A)
        model.CSA(tf_A,ic) = -1;
    end
    tf_B = strcmpi(prot_cx{2},model.Regulators);
    if any(tf_B)
        model.CSA(tf_B,ic) = -1;
    end
    model.CSA(clxind(ic),ic) = 1;
end
model.CSAind = algind;

%---------------------------------------------------------------
n_cmplx = length(find(cplx_ind));
% n_notcmplx = nprot_only +...
%              length(find(rnape_ind)) +...
%              length(find(cm_ind)) +...
%              length(find(metab_ind));
RegCMPLX = nprot_only+1:nprot_only+n_cmplx;
% model.RegCMPLX = n_notcmplx+1:n_notcmplx+n_cmplx;
% model.RegCMPLX_A = zeros(n_cmplx,1);
% model.RegCMPLX_B = zeros(n_cmplx,1);
model.CS = sparse(length(model.Regulators),n_cmplx);
for icplx = 1:n_cmplx
    prot_cx = strsplit(model.Regulators{RegCMPLX(icplx)},'-');
    %Protein/Metabolite Index for Complex (A in A-B complex)
    if ~isempty(prot_cx{1})
        tf_cx1 = strcmpi(prot_cx{1},model.Regulators);
%         model.RegCMPLX_A(icplx) = find(tf_cx1);
        model.CS(tf_cx1,icplx) = -1;
    end
    %Inducer Index for Complex (B in A-B complex)
    if ~isempty(prot_cx{2})
        tf_cx2 = strcmpi(prot_cx{2},model.Regulators);
%         model.RegCMPLX_B(icplx) = find(tf_cx2);
        model.CS(tf_cx2,icplx) = -1;
    end
    model.CS(RegCMPLX(icplx),icplx) = 1;
end
%-------------------------------------------------------------------
nprot_only = nprot_only+n_cmplx;
%Regualtor indices common to both proteins and metabolites 
model.PMind_R = nprot_only + 1:nprot_only + length(find(cm_ind));
model.PMind_M = zeros(1,length(model.PMind_R));
cm_reg = model.Regulators(model.PMind_R);
for icm = 1:length(cm_reg)
    tf_cm = strcmp(cm_reg{icm},model.Metabolites);
    if any(tf_cm)
       model.PMind_M(icm) = find(tf_cm);
    end
end
%Separate phosphorylated proteins from non-phosphorylated proteins
model.PHOSind_M = find(~cellfun('isempty',regexp(model.Metabolites,'^(?:P_)\w')));
model.PHOSind_R = zeros(1, length(model.PHOSind_M));
phos_m = model.Metabolites(model.PHOSind_M);
for iphos = 1:length(phos_m)
    tf_phos = strcmp(phos_m{iphos},model.Regulators);
    if any(tf_phos)
        model.PHOSind_R(iphos) = find(tf_phos);
    end
end
%PTS metabolites
pts_metab = {'pep[c]','pyr[c]','g6p[c]','lac[c]','gl[c]','gal[c]',...
             'glc[e]','lac[e]','gl[e]','gal[c]'};
model.PTSind_M = zeros(length(pts_metab),1);
for ipts = 1:length(pts_metab)
    tf_pts = strcmpi(model.Metabolites,pts_metab{ipts});
    if any(tf_pts)
        model.PTSind_M(ipts) = find(tf_pts);
    end
end
if any(model.PTSind_M)
    model.PTSind_M = model.PTSind_M(model.PTSind_M~=0);
end
model.PTSind_R = zeros(length(pts_metab),1);
for ipts = 1:length(pts_metab)
    tf_pts = strcmpi(model.Regulators,pts_metab{ipts});
    if any(tf_pts)
        model.PTSind_R(ipts) = find(tf_pts);
    end
end
if any( model.PTSind_R )
    model.PTSind_R = model.PTSind_R(model.PTSind_R~=0);
end

%All complexes - Reg and No Reg
%Other Complexes that are part of regulation but do not regulate 
%Eg. LacI-Allolactose in the LAC operon system or cAMP-CRP
cplx = model.Regulators(~cellfun('isempty',strfind(model.Regulators,'-')));
n_clx = length(cplx);
model.CLX = zeros(n_clx,1);
model.CLX_A = zeros(n_clx,1);
model.CLX_B = zeros(n_clx,1);
model.CLX_M = zeros(n_clx,1);
for ic = 1:n_clx
    tf_rc = strcmpi(cplx{icpx},model.Regulators);
    tf_mc = strcmpi(cplx{icpx},model.Metabolites);
    if any(tf_mc)
        %Index of noRegCLX in model.Metabolites
        model.CLX_M(icpx) = find(tf_mc);        
    end
    if any(tf_rc)
        model.CLX(ic) = find(tf_rc);
        model.Kb(tf_rc) = 1e-4;
        model.Kub(tf_rc) = 1e-4;
        cterms = strsplit(cplx{ic},'-');        
        if ~isempty(cterms) && length(cterms) > 1
            if ~isempty(cterms{1})
                %Index for A in complex A-B
                model.CLX_A(ic) = find(strcmpi(cterms{1},model.Regulators));
                if any(strcmpi(cterms{1},model.Regulators))
                    %Active Form of A in A-B                    
%                     model.CLXact_A(icpx) =...
%                     find(strcmpi(sprintf('%s*',cterms{1}),model.Regulators));                
                end                  
            end
            if ~isempty(cterms{2})
                %Index for B in complex A-B
                model.CLX_B(ic) = find(strcmpi(cterms{2},model.Regulators));
            end
        end
    end
end
noregcplx = model.Metabolites(~cellfun('isempty',strfind(model.Metabolites,'-')));
n_noregclx = length(noregcplx);
model.noRegCLX = zeros(n_noregclx,1);
model.noRegCLX_M = zeros(n_noregclx,1);
model.noRegCLX_A = zeros(n_noregclx,1);
model.noRegCLX_B = zeros(n_noregclx,1);
model.noRegCLXact_A = zeros(n_noregclx,1);
% model.noRegCLXact_M = zeros(n_noregclx,1);
for icpx = 1:n_noregclx    
    tf_rc = strcmpi(noregcplx{icpx},model.Regulators);
    tf_mc = strcmpi(noregcplx{icpx},model.Metabolites);
    if any(tf_mc)
        %Index of noRegCLX in model.Metabolites
        model.noRegCLX_M(icpx) = find(tf_mc);        
    end
    if any(tf_rc)
        model.noRegCLX(icpx) = find(tf_rc);
        model.Kb(tf_rc) = 1e-4;
        model.Kub(tf_rc) = 1e-4;
        cterms = strsplit(noregcplx{icpx},'-');        
        if ~isempty(cterms) && length(cterms) > 1
            if ~isempty(cterms{1})
                %Index for A in complex A-B
                model.noRegCLX_A(icpx) = find(strcmpi(cterms{1},model.Regulators));
                if any(strcmpi(cterms{1},model.Regulators))
                    %Active Form of A in A-B                    
                    model.noRegCLXact_A(icpx) =...
                    find(strcmpi(sprintf('%s*',cterms{1}),model.Regulators));                
                end                  
            end
            if ~isempty(cterms{2})
                %Index for B in complex A-B
                model.noRegCLX_B(icpx) = find(strcmpi(cterms{2},model.Regulators));
            end
        end
    end
end

model.nt_gene = length(model.Gene);
model.nt_prot = length(model.Enzyme);
model.nr_cmplx = n_cmplx;
model.nt_metab = length(model.Metabolites);
model.next_metab = FBAmodel.next_metab;
model.nint_metab = nt_metab-model.next_metab;

%Convert matrices to a single parameter vector
model.allpar = parameter_vector(model,length(model.Gene));
model = rmfield(model,{'Coefficient','brate'}); 
model.pmeter = pmeter;

function [ruleterm,tterms] =...
        form_rule(ruleterm,rule,irule,call,tterms,pinducer,minducer,Kb,Kub)
    if nargin < 9
        Kub = {};
    end
    if nargin < 8
        Kb = {};
    end    
    if nargin < 7
        minducer = {};
    end
    if nargin < 6
        pinducer = {};
    end
    opbr = strfind(rule,'(');
    clbr = strfind(rule,')');
    nopbr = length(opbr);    
    newlogic = cell(nopbr-1,1);
    ibr = 1;    
    while irule <= irule+nopbr && ibr <= nopbr
        coeff_terms = rule(opbr(ibr):clbr(ibr));        
        orterms = strsplit(rule(opbr(ibr)+1:clbr(ibr)-1),'|');
        andterms = strsplit(rule(opbr(ibr)+1:clbr(ibr)-1),'&');
        andpos = strfind(rule(opbr(ibr)+1:clbr(ibr)-1),'&');        
        if ~isemptyr(orterms) && isempty(andpos)
            relcoeff = coeff(tterms+1:tterms+length(orterms));
            tterms = tterms+length(orterms);
            if ~assign_coeff || length(relcoeff) ~= length(orterms)
                [relcoeff,coeff_flag] = coeffinit(coeff_terms,defparval,coeff_flag);            
            end
            [model,ruleterm] =...
            assignterm(model,ruleterm,irule+1,orterms,'|',call,relcoeff,...
                       pinducer,minducer,Kb,Kub);
        elseif ~isemptyr(andterms) 
            relcoeff = coeff(tterms+1:tterms+length(andterms));
            tterms = tterms+length(andterms);
            if ~assign_coeff || length(relcoeff) ~= length(andterms)
                [relcoeff,coeff_flag] = coeffinit(coeff_terms,defparval,coeff_flag);
            end
            [model,ruleterm] =...
            assignterm(model,ruleterm,irule+1,andterms,'&',call,relcoeff,...
                       pinducer,minducer,Kb,Kub);
        end
        if ibr <= nopbr-1
            newlogic{ibr} = rule(clbr(ibr)+1:opbr(ibr+1)-1); 
            %irule = irule+1;
        end
        if ibr < nopbr && ~isempty(newlogic{ibr})
            ruleterm{irule+1} = [ruleterm{irule+1},newlogic{ibr}];
            %irule = irule+1;
        end        
        ibr = ibr+1;        
        irule = irule+1;
    end        
end
    
function [model,rule] =...
        assignterm(model,rule,ibr,terms,sym,call,relcoeff,pinducer,minducer,kb,kub)
    jterm = 1;
    kprot = length(model.Enzyme)+1;
    jmetab = length(model.Metabolites)+1;
    jreg = length(model.Regulators)+1;
    if ~isempty(pinducer) || ~isempty(minducer)        
        p_terms = strsplit(strtrim(strrep(pinducer,'"','')),',');
        m_terms = strsplit(strtrim(strrep(minducer,'"','')),',');
        try
            inducers = [p_terms m_terms];
        catch
            inducers = [p_terms;m_terms];
        end   
        inducers = inducers(~cellfun('isempty',inducers));
    end
    relKb = zeros(length(terms),1);
    relKub = zeros(length(terms),1);
    while jterm <= length(terms)
        [match] = regexp(terms{jterm},'(\S+.?)(\W+.?)','tokens');
%         [match] = regexp(terms{jterm},'(\w+.?)\[(\W+.?)\]','tokens');   
        if isempty(match)
            fprintf('\nMissing regulator sign\nCheck Data File\n Gene:%d',igene);
%             term_flag = 1;
            return
        end
        switch match{1}{2}
            case '+'
                stoich = 1;
            case '-'
                stoich = -1;
            case '+/-'
                stoich = 2;
        end
        %Building Rules        
        if call == 1
            if jterm == 1 && length(rule) >= ibr
                rule{ibr} = [rule{ibr},match{1}{1}];
            elseif jterm == 1
                rule = [rule;match{1}{1}];
            elseif jterm > 1
                rule{ibr} = [rule{ibr},sprintf('%s%s',sym,match{1}{1})];
            end
            %Proteins Only
            tfp = strcmpi(match{1}{1},model.Enzyme);
            if ~any(tfp) && (isempty(pinducer)&&isempty(minducer))
                model.Enzyme{kprot} = match{1}{1};
                %model.Prot(kprot) = 0;
                model.SSprot(kprot) = 0;
                kprot = kprot+1;                
            elseif ~any(tfp) && (~isempty(pinducer)||~isempty(minducer))
                %Protein-Inducer/Operator Complex                
                prot_cmp = strsplit(match{1}{1},'-');                    
                %Protein part
                tfp_c = strcmpi(prot_cmp{1},model.Enzyme);
                if ~any(tfp_c)
                    model.Enzyme{kprot} = prot_cmp{1};                
                    model.SSprot(kprot) = 0;
                    kprot = kprot+1;
                end
                %Assigning Protein as Regulators
                [model,jreg] = assgIndReg_(prot_cmp{1},model,jreg);
                if length(prot_cmp)>1 %If inducer complex is present
                    %Inducer part
                    [model,kprot,jmetab] =...
                    inducer_(model,prot_cmp{2},pinducer,minducer,kprot,jmetab);
                    %Assigning Inducer as Regulator
                    [model,jreg] = assgIndReg_(prot_cmp{2},model,jreg);
                    %Assigning Kbing and Kunbinding
                    [relKb(jterm),relKub(jterm)] = inducerK(kb,kub,prot_cmp{2},inducers);
                end
                %Assigning Complex as Regulator
                %Done outside the if-else call                
            end              
        elseif call == 2
            if iscell(rule)
                rule = changeGeneRules(rule,match{1}{1},sym);
            else
                rule = changeGeneRules({rule},match{1}{1},sym);
            end
            %Metabolites Only
            tfm = strcmpi(match{1}{1},model.Metabolites);
            if ~any(tfm) && (isempty(pinducer)&&isempty(minducer))
                model.Metabolites{jmetab} = match{1}{1};
                model.Metab(jmetab) = 0;
                model.maxMetab(jmetab) = 0;
                jmetab = jmetab+1;
            elseif ~any(tfm) && (~isempty(pinducer)||~isempty(minducer))
                %Metabolite-Inducer Complex
                met_cmp = strsplit(match{1}{1},'-');
                %Metabolite part
                tfm_c = strcmpi(met_cmp{1},model.Metabolites);
                if ~any(tfm_c)
                    model.Enzyme{kprot} = met_cmp{1};                
                    model.SSprot(kprot) = 0;
                    jmetab = jmetab+1;
                end
                %Assigning Metabolites as Regulators
                [model,jreg] = assgIndReg_(met_cmp{1},model,jreg);
                if length(met_cmp)>1
                    %Inducer Part
                    [model,kprot,jmetab] =...
                    inducer_(model,met_cmp{2},pinducer,minducer,kprot,jmetab);
                    %Assigning Inducer as Regulator
                    [model,jreg] = assgIndReg_(met_cmp{2},model,jreg);
                    %Assigning Kbing and Kunbinding
                    [relKb(jterm),relKub(jterm)] = inducerK(kb,kub,prot_cmp{2},inducers);
                end
                %Assigning Complex as Regulator
                %Done outside the if-else call  
            end
        end

        tf_all = strcmpi(match{1}{1},model.Regulators);           
        if any(tf_all)
            model.RS(igene,tf_all) = stoich;    
            try
                model.Coefficient(igene,tf_all) = str2double(relcoeff{jterm}{1});
            catch
                model.Coefficient(igene,tf_all) = str2double(relcoeff{jterm});
            end       
            model.Kb(tf_all) = relKb(jterm);%default Kbind values
            model.Kub(tf_all) = relKub(jterm);%default Kubind values
        else            
            model.Regulators{jreg} = match{1}{1};
            model.RS(igene,jreg) = stoich;                               
            model.Coefficient(igene,jreg) = str2double(relcoeff{jterm}{1});
            model.trate(jreg,igene) = 0;
            model.Kb(jreg) = relKb(jterm);%default Kbind values
            model.Kub(jreg) = relKub(jterm);%default Kubind values
            model.relReg(jreg) = 0;
            model.SSreg(jreg) = 0;            
            jreg = jreg + 1;
        end    
        jterm = jterm + 1;
    end
end

function [model,kprot,jmetab] =...
        inducer_(model,inducer,pinducer,minducer,kprot,jmetab)    
    %Assign Inducer to Protein or Metabolite 
    if ~isempty(pinducer)%Protein Inducer
        if ~isempty(inducer)
            tfp_c = strcmpi(inducer,model.Enzyme);
            if ~any(tfp_c)
                model.Enzyme{kprot} = inducer;                
                model.SSprot(kprot) = 0;
                kprot = kprot+1;                
            end
        end
    end
    if ~isempty(minducer)%Metabolite Inducer
        if ~isempty(inducer)
            tfm_c = strcmpi(inducer,model.Metabolites);
            if ~any(tfm_c)
                model.Metabolites{jmetab} = inducer;
                model.Metab(jmetab) = 0;
                model.maxMetab(jmetab) = 0;
                jmetab = jmetab+1;
            end
        end
    end 
end
function [model,jreg] = assgIndReg_(species,model,jreg)
    %Assign Species "species" as Regulator
    tf_all = strcmpi(species,model.Regulators);
    if ~any(tf_all)
        model.Regulators{jreg} = species;
        model.RS(igene,jreg) = 0;                                   
        model.Coefficient(igene,jreg) = 0;
        model.trate(jreg,igene) = 0;
        model.relReg(jreg) = 0;
        model.SSreg(jreg) = 0;
        model.Kb(jreg) = 0;
        model.Kub(jreg) = 0;
        jreg = jreg + 1;
    end
end
function [relKb,relKub] = inducerK(kb,kub,inducer,inducer_list)
    tf_pind = strcmpi(inducer,inducer_list);
    if any(tf_pind)
        relKb = str2double(kb{tf_pind});
        relKub = str2double(kub{tf_pind});
    end
end
%==========================================================================
%Regulatory Protein Concentration
%trnmodel.ProtConc = ones(length(trnmodel.Protein),1);
%Converting # molecules to Concentration
%Parameters
%Cell Volume = 1E-15 L
%Na = 6.023E+23 molecules/mole
%Concentration = M mole/L
%Order of magnitude of Concentrations = 1E-09 (nM range)
%trnmodel.ProtConc = trnmodel.ProtConc/(1E-15*6.023E+23);
%==========================================================================
%Translation Rates for Gene-based Proteins
%See addgene.m for changes 
%trnmodel.trate is populated there
%November 2013
end
