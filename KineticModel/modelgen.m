function [model_data,parameter,variable,nt_rxn,nt_metab,bmrxn] = modelgen(rxfname)
%Generate Metabolic Network Model from Stoichiometric Reactions in the
%Metabolic Network

%Shyam 2014
fprintf('Creating MATLAB model from model file...\n');
fileid = fopen(rxfname);
if fileid == -1
    fprintf('File %s cannot be opened.', rxfname);
    model_data = struct([]);
    parameter = struct([]);
    variable = struct([]);
    return;
end

C = textscan(fileid, '%s%s%s%s%f%f%f%f%f%f%f%f%f%s%s%s%s%s%s%f%f',...
                     'Delimiter', '\t',...
                     'TreatAsEmpty', {'None'},...
                     'HeaderLines', 1);
fclose(fileid);

model_data = struct();
model_data.nt_rxn = length(C{3}(~cellfun('isempty',C{3})));

nt_rxn = model_data.nt_rxn;
%enzs/Reaction Information 
model.enzName = C{1}(~cellfun('isempty',C{1}));%enzs Name
%model.RxnMech = C{12};
model.S = sparse(0,nt_rxn);
model.K = sparse(0,nt_rxn);
model.delG = zeros(nt_rxn,1);
model.kcat_fwd = C{9}(1:nt_rxn);   %Reaction fwd kcat 
model.kcat_bkw = C{10}(1:nt_rxn);   %reaction bkw kcat
model.Vmax = C{12}(1:nt_rxn);   %Vmax
model.reversible = zeros(nt_rxn,1);
model.Vss = zeros(nt_rxn,1);
model.mets = {};
Keq = zeros(nt_rxn,1);
model.EC = C{11}(1:nt_rxn);     %Etotal Concentration
model.Klb = sparse(0,nt_rxn);
model.Kub = sparse(0,nt_rxn);
model.delGlb = zeros(nt_rxn,1);
model.delGub = zeros(nt_rxn,1);

imetab = 1;
for irxn = 1:nt_rxn
    ipt = 0; 
    %Vsteady state
    if ~isempty(C{13}(irxn)) && ~isnan(C{13}(irxn))
        model.Vss(irxn) = C{13}(irxn);
    else
        %run FBA after building model
        %runFBA = true;
        model.Vss(irxn) = 1;
    end   
    
    %delG bounds
    if ~isempty(C{6}(irxn)) && ~isnan(C{6}(irxn))
        model.delGlb(irxn) = C{6}(irxn);
    else
        model.delGlb(irxn) = 0;    
    end
    
    if ~isempty(C{7}(irxn)) && ~isnan(C{7}(irxn))
        model.delGub(irxn) = C{7}(irxn);
    else
        model.delGub(irxn) = 0;
    end
    
    %Keq if available else use delG
    %use delG to obtain Keq
    %delG = -RTlnKeq
    if ~isempty(C{8}(irxn)) && ~isnan(C{8}(irxn))
        
        Keq(irxn) = C{8}(irxn);
        if ~isempty(C{5}(irxn)) && ~isnan(C{5}(irxn))
            model.delG(irxn) = C{5}(irxn);
        else
            model.delG(irxn) = -0.008314*298*log(Keq(irxn));
        end
        
    elseif ~isempty(C{5}(irxn)) && ~isnan(C{5}(irxn))
        
        model.delG(irxn) = C{5}(irxn);
        if ~isempty(C{8}(irxn)) && ~isnan(C{8}(irxn))
            Keq(irxn) = C{8}(irxn);
        else
            Keq(irxn) = exp(-C{5}(irxn)/(0.008314*298.15));
        end

    end
    
    %Building S and K matrices
    rxnstring = C{3}{irxn};  
    %separate terms into a vector
    [par,Klb,Kub] = extract_par(C{15}{irxn});    
    
    if ~isempty(strfind(rxnstring, '<==>'))
		eqsym = '<==>';
		reverse = 1;
    elseif ~isempty(strfind(rxnstring,'<=>'))
        eqsym = '<=>';
		reverse = 1;
    elseif ~isempty(strfind(rxnstring, '--->'))
		eqsym = '--->';
		reverse = 0;
    elseif ~isempty(strfind(rxnstring,'-->'))
        eqsym = '-->';
        reverse = 0;
    elseif ~isempty(strfind(rxnstring, '->'))
		eqsym = '->';
		reverse = 0;
    elseif ~isempty(strfind(rxnstring, '='))
		eqsym = '=';
		reverse = 1;
    end    

    k = strfind(rxnstring, ':');
    if ~isempty(k)
        compartment = rxnstring(1:(k - 2));
        rxnstring = rxnstring((k + 1):end);
    else
        compartment = '';
    end
    
    %Reaction LHS
    k = strfind(rxnstring, eqsym);
    lhs = rxnstring(1:k - 1);
    if ~isempty(lhs)
        terms = strtrim(textscan(lhs, '%s', 'Delimiter', '+'));
        s = regexp(terms{1}, '[(]?([0-9.]+)[)]? ([A-Za-z0-9_\-\[\]]+)', 'tokens');  
        %Assign default parameters here if isempty(par) == 1
        if isempty(par)
            %call function to assign default parameters for lhs
            par = defparval(length(s));
        elseif any(par == 0)
            par(par == 0) = defparval(length(find(par==0)));
        end

        for iterm = 1:length(s)
            if ~isempty(s{iterm})
                stoich = str2double(s{iterm}{1}{1});
                if ~isempty(compartment)
                    metab = [s{iterm}{1}{2} compartment];
                else
                    metab = s{iterm}{1}{2};
                end
            else
                stoich = 1;
                if ~isempty(compartment)
                    metab = [terms{1}{iterm} compartment];
                else
                    metab = terms{1}{iterm};
                end
            end   
            
            %Required for compartmentalized models
            tf = strcmpi(metab, model.mets);
            if any(tf)
                model.S(tf, irxn) = -stoich;                
                if ~isempty(par)
                    model.K(tf,irxn) = par(ipt+iterm);
                    model.Klb(tf,irxn) = 0;
                    model.Kub(tf,irxn) = 0;
                end
                if ~isempty(Klb)
                    model.Klb(tf,irxn) = Klb(ipt+iterm);
                end
                if ~isempty(Kub)
                    model.Kub(tf,irxn) = Kub(ipt+iterm);
                end
            else
                model.S(imetab, irxn) = -stoich;                
                model.mets{imetab} = metab;
                if ~isempty(par)
                    model.K(imetab,irxn) = par(ipt+iterm);
                    model.Klb(imetab,irxn) = 0;
                    model.Kub(imetab,irxn) = 0;
                end
                if ~isempty(Klb)
                    model.Klb(imetab,irxn) = Klb(ipt+iterm);
                end
                if ~isempty(Kub)
                    model.Kub(imetab,irxn) = Kub(ipt+iterm);
                end
                imetab = imetab + 1;
            end
        end
        ipt = ipt + iterm;
    end
    
    %rwaction RHS
    rhs = rxnstring((k + length(eqsym)):end);
    if ~isempty(rhs)
        terms = strtrim(textscan(rhs, '%s', 'Delimiter', '+'));
        s = regexp(terms{1},'[(]?([0-9.]+)[)]? ([A-Za-z0-9_\-\[\]]+)','tokens');
        if length(par) < ipt + length(s) 
            %-> assign default parameters for rhs
            par = defparval(length(s),par);
        end
        for iterm = 1:length(s)
            if ~isempty(s{iterm})
                stoich = str2double(s{iterm}{1}{1});
                if ~isempty(compartment)
                    metab = [s{iterm}{1}{2} compartment];
                else
                    metab = s{iterm}{1}{2};
                end
            else
                stoich = 1;
                if ~isempty(compartment)
                    metab = [terms{1}{iterm} compartment];
                else
                    metab = terms{1}{iterm};
                end
            end            
            tf = strcmpi(metab, model.mets);
            if any(tf)
                model.S(tf, irxn) = stoich;                
                if ~isempty(par)
                    model.K(tf,irxn) = par(ipt+iterm);
                    model.Klb(tf,irxn) = 0;
                    model.Kub(tf,irxn) = 0;
                end
                if ~isempty(Klb)
                    model.Klb(tf,irxn) = Klb(ipt+iterm);
                end
                if ~isempty(Kub)
                    model.Kub(tf,irxn) = Kub(ipt+iterm);
                end
            else
                model.S(imetab, irxn) = stoich;                
                model.mets{imetab} = metab;
                if ~isempty(par)
                    model.K(imetab,irxn) = par(ipt+iterm);
                    model.Klb(imetab,irxn) = 0;
                    model.Kub(imetab,irxn) = 0;
                end
                if ~isempty(Klb)
                    model.Klb(imetab,irxn) = Klb(ipt+iterm);
                end
                if ~isempty(Kub)
                    model.Kub(imetab,irxn) = Kub(ipt+iterm);
                end
                imetab = imetab + 1;
            end
        end
        ipt = ipt + iterm;
    end 
    model.reversible(irxn) = reverse;           
end
model.mets = model.mets';
model.nt_metab = size(model.S, 1);
nt_metab = model.nt_metab;

% Selecting activators/inhibitors and building SI and KI
model.SI = sparse(nt_metab,nt_rxn);%Regulatory S matrix
model.SItype = sparse(nt_metab,nt_rxn);%Regulation type
% model.KI = sparse(nt_metab,nt_rxn);
model.KIact = sparse(nt_metab,nt_rxn);
model.KIihb = sparse(nt_metab,nt_rxn);
for irxn = 1:nt_rxn
    [par,Klb,Kub] = extract_par(C{18}{irxn});%Acquire parameters as vectors
    ireg = 0;       
    actstring = strtrim(strrep(C{16}{irxn},'"',''));%Activators    
    if ~isempty(actstring)         
        [model] = ident_regulator(model,actstring,1,par,Klb,Kub);%Activators        
    end    
    inhstring = strtrim(strrep(C{17}{irxn},'"',''));%Inhibitors
    if ~isempty(inhstring)  
        %->Assign default parameters for inhibitors if par = [] or 
        %if length(par) < length(activators) + length(inhibitors)
        %par = defparval(nterms,par)
        [model] = ident_regulator(model,inhstring,-1,par,Klb,Kub);%Inhibitors        
    end
end

% Separate External & Internal mets
% exter_mind = ~cellfun('isempty',regexp(model.mets,'\w(?:xt)$'));
newmodel = separate_cex(model);
model_data.mets = newmodel.mets;                   

% Identify Reaction Indices
model_data.rxns = model.enzName;
[Vind,VFex,Vex,bmrxn] = fluxIndex(model_data,nt_rxn,newmodel.S);

% Append appropriate rows/columns corresponding to enzss/fluxes
nenz = length(model.enzName);
other_ind = setdiff(1:nenz,[ToColumnVector(Vind)...                            
                            ToColumnVector(Vex)...
                            ToColumnVector(VFex)...
                            bmrxn]);
[mS,~] = size(newmodel.S);
model_data.enzs = [model.enzName(setdiff(1:nenz,other_ind),1);...
                     model.enzName(other_ind,1)];
model_data.rxns = model_data.enzs;    
newmodel.S = [newmodel.S,sparse(mS,length(other_ind))];
newmodel.SI = [newmodel.SI,sparse(mS,length(other_ind))];
newmodel.K = [newmodel.K,sparse(mS,length(other_ind))];
newmodel.Klb = [newmodel.Klb,sparse(mS,length(other_ind))];
newmodel.Kub = [newmodel.Kub,sparse(mS,length(other_ind))];
newmodel.KIact = [newmodel.KIact,sparse(mS,length(other_ind))];
newmodel.KIihb = [newmodel.KIihb,sparse(mS,length(other_ind))];
newVss = [model.Vss;zeros(length(other_ind),1)];
delSGr = [model.delG;zeros(length(other_ind),1)];
delGlb = [model.delGlb;zeros(length(other_ind),1)];
delGub = [model.delGub;zeros(length(other_ind),1)];
newKeq = [Keq;zeros(length(other_ind),1)];
newkcfwd = [model.kcat_fwd;zeros(length(other_ind),1)];
newkcbkw = [model.kcat_bkw;zeros(length(other_ind),1)];
newreverse = [model.reversible;zeros(length(other_ind),1)];

% compensated metabolites
model_data.CMPS = sparse(length(model_data.mets),length(model_data.rxns));
% compensated species
for irxn = 1:nt_rxn
    rxnstring = C{3}{irxn};
    k = strfind(rxnstring, ':');
    if ~isempty(k)
        compartment = rxnstring(1:(k - 2));        
    else
        compartment = '';
    end
    cmp_sp = C{4}{irxn};
    if ~isempty(cmp_sp)
        cmp_sp = strtrim(strrep(cmp_sp,'"',''));
        cmp_sp = strsplit(cmp_sp,{' ',','},'CollapseDelimiters',true);
        for isp = 1:length(cmp_sp)
            if ~isempty(compartment)
                cmp_sp{isp} = [cmp_sp{isp} compartment];
                tf = strcmpi(model_data.mets,cmp_sp{isp});
                if any(tf)
                    model_data.CMPS(tf,irxn) = 1;
                end
            else
                tf = strcmpi(model_data.mets,cmp_sp{isp});
                if any(tf)
                    model_data.CMPS(tf,irxn) = 1;
                end
            end
        end
    end
end

% New Indices
[Vind,VFex,Vex,bmrxn] = fluxIndex(model_data,nt_rxn,newmodel.S);
% [Vind,Vuptake,VFup,VFex,Vex,bmrxn,Vup,Vdn] = fluxIndex(model_data,nt_rxn,newS);
model_data.Vind = Vind;
model_data.Vex = Vex;
model_data.VFex = VFex;
model_data.bmrxn = bmrxn;
% Identify activated reactions
[~,allactrxn] = find(newmodel.SI(:,1:nt_rxn)>0);
Vact_ind = unique(allactrxn);
% Identify Inhibited reactions
[~,allihbrxn] = find(newmodel.SI(:,1:nt_rxn)<0);
Vihb_ind = unique(allihbrxn);
% K and KI are written as vectors as opposed to matrices
parameter = struct();
model_data.S = newmodel.S;
% model_data.CMPS = newCMPS;
model_data.SI = newmodel.SI;
model_data.Vss = newVss;
model_data.delSGr = delSGr;
model_data.delGlb = delGlb;
model_data.delGub = delGub;
model_data.Keq = newKeq;
model_data.rev = newreverse;

model_data.Vact_ind = Vact_ind;
model_data.Vihb_ind = Vihb_ind;
model_data.bmrxn = bmrxn;

model_data.n_rxn = length(model_data.Vind);
model_data.nt_metab = nt_metab;
model_data.next_metab = newmodel.next_metab;
model_data.nint_metab = newmodel.nint_metab;

parameter.K = newmodel.K;
parameter.Klb = newmodel.Klb;
parameter.Kub = newmodel.Kub;
parameter.KIact = newmodel.KIact;
parameter.KIihb = newmodel.KIihb;
parameter.Vmax = model.Vmax;
parameter.kcat_fwd = newkcfwd;
parameter.kcat_bkw = newkcbkw;

%Separate phosphorylated proteins from non-phosphorylated proteins
%Separate PTS mets from other mets
%At least get indices corresponding to PTS mets
%PEP,PYR,G6P,GLCxt,ACxt,RIBxt,GLXxt, etc
pts_metab = {'pep[c]','pyr[c]','g6p[c]','lac[c]','gl[c]','gal[c]',...
             'glc[e]','lac[e]','gl[e]','gal[e]'};
pts_ind = zeros(length(pts_metab),1);                                                                                                                                                                                                                                                                                                                                                                                          
for ipts = 1:length(pts_metab)
    tf_pts = strcmpi(model_data.mets,pts_metab{ipts});
    if any(tf_pts)
        pts_ind(ipts) = find(tf_pts);
    end
end
if any(pts_ind)
    pts_ind = pts_ind(pts_ind~=0);
end
model_data.PTSind = pts_ind;

% Reading in the concentration for each metabolite
metabname = C{19}(~cellfun('isempty',C{19}));%Metabolite Name
concLow = C{20}(1:length(metabname));%Metabolite concentration
concHigh = C{21}(1:length(metabname));
model_data.MClow = zeros(nt_metab,1);
model_data.MChigh = zeros(nt_metab,1);
for imc = 1:length(model_data.mets)
    mtf = strcmpi(metabname,model_data.mets{imc});
    if any(mtf) 
        if ~isnan(concLow(mtf))
            model_data.MClow(imc) = concLow(mtf);        
        else
            model_data.MClow(imc) = 0.001;%Randon Concentration for id purposes            
        end
        if ~isnan(concHigh(mtf))
            model_data.MChigh(imc) = concHigh(mtf);
        else
            model_data.MChigh(imc) = 5;
        end
    else%No concentrations are specified in file
        model_data.MClow(imc) = 0.001;
        model_data.MChigh(imc) = 5;
    end
end
model_data.b = zeros(nt_metab,1);
model_data.c = sparse(1,bmrxn,1,1,nt_rxn)';
model_data.vl = zeros(nt_rxn,1);
model_data.vl(model_data.vl==0) = -500;
model_data.vl(bmrxn) = 0;
model_data.vu = zeros(nt_rxn,1);
model_data.vu(model_data.vu==0) = 500;

%if runFBA
    %build FBA matrices
    %set flux bounds
    %set uptake rates
    %run FBA
    %obtain ss fluxes
%end

variable = struct();
% variable.MC = model.MC;
variable.EC = model.EC;

%Add Molecular weight individually
model_data.MolWt = zeros(nt_metab,1);
model_data.MolWt(strcmpi('g3p[c]',model_data.mets)) = 172.074;
model_data.MolWt(strcmpi('pyr[c]',model_data.mets)) = 88.06;
model_data.MolWt(strcmpi('pep[c]',model_data.mets)) = 168.042;
model_data.MolWt(strcmpi('f6p[c]',model_data.mets)) = 259.81;
model_data.MolWt(strcmpi('g6p[c]',model_data.mets)) = 260.136;
model_data.MolWt(strcmpi('B[c]',model_data.mets)) = 200;

fprintf('Model generation complete\n\n');
% nested functions
function [model] =...
ident_regulator(model,reg_string,reg_stoich,par,Klb,Kub)%pass par as argument
    %KI to added 
    compos = strfind(reg_string,','); 
    regterms = cell(length(compos)+1,1); 
    if ~isemptyr(compos)
        for icompos = 1:length(compos) 
            if icompos == 1
                regterms{icompos} = reg_string(1:compos(icompos)-1);
            else
                regterms{icompos} =...
                reg_string(compos(icompos-1)+1:compos(icompos)-1); 
            end
        end
        regterms{icompos+1} = reg_string(compos(icompos)+1:end);
    else
        regterms{1} = reg_string(1:end);
    end    
    nterms = length(regterms);
    
    %->Assign default parameters 
    if isempty(par) 
        par = defparval(nterms);
    elseif length(par) <= ireg
        par = defparval(nterms,par);
    elseif any(par==0)
        par(par == 0) = defparval(length(find(par==0)));
    end

    iregterm = 1;    
    while iregterm <= nterms
        [mech,mechx] = regexp(regterms{iregterm},'(\w+.?)\((\w+.?)\)+','tokens','split');
        mechx = mechx(~cellfun('isempty',mechx));
        if ~isempty(mechx)
            newterms = regexp(mechx{1},'(\w+.?)+','tokens');
            mech = [mech,newterms];
        end
        if iscell(mech{1})
            mech = mech{1};
        end

        metabindx = strcmpi(mech{1},model.mets);            
        if any(metabindx)
            model.SI(metabindx,irxn) = reg_stoich;
            if ~isempty(par)
%                     model.KI(metabindx,irxn) = par(ireg+iregterm);
                if reg_stoich > 0
                    model.KIact(metabindx,irxn) = par(ireg+iregterm);
                elseif reg_stoich < 0
                    model.KIihb(metabindx,irxn) = par(ireg+iregterm);
                end
                model.Klb(metabindx,irxn) = 0;
                model.Kub(metabindx,irxn) = 0;
            end
            if ~isempty(Klb)
                model.Klb(metabindx,irxn) = Klb(ireg+iregterm);
            end
            if ~isempty(Kub)
                model.Kub(metabindx,irxn) = Kub(ireg+iregterm);
            end
            if length(mech) < 2%no mechanism specified
                [model] = reg_type('O',[find(metabindx);irxn],model);                        
            else
                [model] = reg_type(mech{2},[find(metabindx);irxn],model);
            end
        else
            model.SI(imetab,irxn) = reg_stoich;
            model.S(imetab,:) = sparse(1,size(model.S,2));
            model.K(imetab,:) = sparse(1,size(model.K,2));
            model.mets{imetab} = mech{1};
            if ~isempty(par)
%                     model.KI(imetab,irxn) = par(ireg+iregterm);
                if reg_stoich > 0
                    model.KIact(imetab,irxn) = par(ireg+iregterm);
                    model.KIihb(imetab,:) = sparse(1,size(model.KIihb,2));
                elseif reg_stoich < 0
                    model.KIihb(imetab,irxn) = par(ireg+iregterm);
                    model.KIact(imetab,:) = sparse(1,size(model.KIact,2));
                end
                model.Klb(imetab,irxn) = 0;
                model.Kub(imetab,irxn) = 0;
            end
            if ~isempty(Klb)
                model.Klb(imetab,irxn) = Klb(ireg+iregterm);
            end
            if ~isempty(Kub)
                model.Kub(imetab,irxn) = Kub(ireg+iregterm);
            end                
            if length(mech) < 2%no mechanism specified
                [model] = reg_type('O',[imetab;irxn],model);                        
            else
                [model] = reg_type(mech{2},[imetab;irxn],model);
            end                    
            imetab = imetab + 1;
        end 
        iregterm = iregterm + 1;
    end  
    ireg = ireg + iregterm - 1;   
end


% function [parameter,lb,ub] = extract_par(par_string)
%     %separate terms into a vector
%     if ~isempty(par_string)
%         par_string = strtrim(strrep(par_string,'"',''));        
%         parm = strsplit(par_string,',');                      
%         parameter = zeros(length(parm),1);
%         lb = zeros(length(parm),1);
%         ub = zeros(length(parm),1);
%         if ~isempty(parm)
%             ipar = 1;
%             while ipar <= length(parm)
%                 parm{ipar} = strtrim(parm{ipar});
%                 bnd_op = strfind(parm{ipar},'[');
%                 bnd_cl = strfind(parm{ipar},']');
%                 if ~isempty(bnd_op) && ~isempty(bnd_cl)
%                     %if there are bounds specified                    
%                     Klub = strsplit(parm{ipar}(bnd_op+1:bnd_cl-1));                        
%                     lb(ipar) = str2double(strtrim(Klub{1}));
%                     ub(ipar) = str2double(strtrim(Klub{2}));
%                 else%No brackets
%                     Kpar = str2double(parm{ipar});
%                     if Kpar ~= 1
%                         parameter(ipar) = Kpar;                            
%                     end                            
%                 end
%                 ipar = ipar + 1;
%             end
%         end
%     else
%         parameter = [];
%         lb = [];
%         ub = [];
%         %Default coefficients are defined within respective modules        
%     end
% end

function [model] = reg_type(mechanism,index,model)
    switch mechanism
        case 'A'
            model.SItype(index(1),index(2)) = 1;                    
        case 'C'
            model.SItype(index(1),index(2)) = 2;
        case 'U'
            model.SItype(index(1),index(2)) = 3;
        case 'N'
            model.SItype(index(1),index(2)) = 4;
        otherwise
            model.SItype(index(1),index(2)) = 5;
    end
end
end
