function [model] = buildmat(fname)
%Build S and Z matrices for Injectivity-based preclusion of
%multistationarity in biological newtorks
%Algorithm based on Feliu and Wiuf 2013 and Feliu and Wiuf 2012
%This function is called from Mathematica 10.0 for symbolic computation of
%determinants
%Shyam 2015
fileid = fopen(fname);
if fileid == -1
    fprintf('File %s cannot be opened.', fname);
    model.S = [];
    model.Z = [];
    return;
end

C = textscan(fileid, '%s%s%s%s%s',...
                     'Delimiter', '\t',...
                     'TreatAsEmpty', {'None'},...
                     'HeaderLines', 1);
fclose(fileid);

nt_rxn = length(C{1}); 
model.S = sparse(0,nt_rxn);
model.Z = sparse(nt_rxn,0);
model.metab = {};
% poseff = cell(nt_rxn);
% negeff = cell(irxn);

imetab = 1;
for irxn = 1:nt_rxn
    nts_rxn = size(model.S,2);
    
    %Effectors to build influence matrix
    poseff_f = effectors(C{2}{irxn});
    negeff_f = effectors(C{3}{irxn});
    poseff_b = effectors(C{4}{irxn});
    negeff_b = effectors(C{5}{irxn});
    
    %Building S and Z matrices
    rxnstring = C{1}{irxn};  
    %separate terms into a vector

    if ~isempty(strfind(rxnstring, '<==>'))
		eqsym = '<==>';
		reverse = 1;
    elseif ~isempty(strfind(rxnstring, '-->'))
		eqsym = '--->';
		reverse = 0;
    elseif ~isempty(strfind(rxnstring,'--->'))
        eqsym = '--->';
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
        %compartment = rxnstring(1:(k - 2));
        rxnstring = rxnstring((k + 1):end);
    else
        %compartment = '';
    end
     
    k = strfind(rxnstring, eqsym);
    lhs = rxnstring(1:k - 1);
    rhs = rxnstring((k + length(eqsym)):end);
    if reverse
        lhs_r = rhs;
        rhs_r = lhs;
    else
        lhs_r = [];
        rhs_r = [];
    end
    if ~isempty(lhs)
        stoich = -1;
        [model,imetab] = rxn_sp(lhs,imetab,irxn,stoich,model,poseff_f,negeff_f);
    end    
    if ~isempty(rhs)
        stoich = 1;
        [model,imetab] = rxn_sp(rhs,imetab,irxn,stoich,model,poseff_f,negeff_f);
    end
    if ~isempty(lhs_r)
        stoich = -1;
        [model,imetab] = rxn_sp(lhs_r,imetab,nts_rxn+1,stoich,model,poseff_b,negeff_b);
    end
    if ~isempty(rhs_r)
        stoich = 1;
        [model,imetab] = rxn_sp(rhs_r,imetab,nts_rxn+1,stoich,model,poseff_b,negeff_b);
    end
%     if ~isempty(lhs_r) && ~isempty(rhs_r)
%         nt_rxn = nt_rxn + 1;
%     end
end
nts_rxn = size(model.S,2);
if size(model.Z,2)<size(model.S,1)
    model.Z = [model.Z zeros(size(model.Z,1),size(model.S,1)-size(model.Z,2))];
end
if size(model.Z,1)<size(model.S,2)
    model.Z = [model.Z;zeros(size(model.S,2)-size(model.Z,1),size(model.Z,2))];
end
% model.Z = sparse(nts_rxn,length(model.metab));
% for irxn = 1:nt_rxn
%     %Effectors to build influence matrix
%     poseff = effectors(C{2}{irxn});
%     negeff = effectors(C{3}{irxn});
%     if ~isempty(poseff)
%         ipeff = 1;
%         while ipeff <= length(poseff)
%             tfpos = strcmp(poseff{ipeff},model.metab);
%             if any(tfpos)
%                 model.Z(irxn,tfpos) = 1;
%             end
%             ipeff = ipeff + 1;
%         end
%     end
%     if ~isempty(negeff)
%         ineff = 1;
%         while ineff <= length(negeff)
%             tfneg = strcmp(negeff{ineff},model.metab);
%             if any(tfneg)
%                 model.Z(irxn,tfneg) = -1;
%             end
%             ineff = ineff + 1;
%         end
%     end        
% end
% tfSt = strcmp('St',model.metab);
% tfS = strcmp('S',model.metab);
% tfX = strcmp('X',model.metab);
% tfE = strcmp('E',model.metab);
% if any(tfSt)&& any(tfS) && any(tfX) && any(tfE)
%     model.S = [model.S(tfS,:);model.S(tfSt,:);model.S(tfE,:);model.S(tfX,:)];
%     model.Z = [model.Z(:,tfS) model.Z(:,tfSt) model.Z(:,tfE) model.Z(:,tfX)];
% end
S = model.S;
Z = model.Z;
model.metab = model.metab';
% model.S = full(model.S);
% model.Z = full(model.Z);
return

function [model,imetab] = rxn_sp(sptrm,imetab,irxn,stmul,model,poseff,negeff)
terms = strtrim(textscan(sptrm, '%s', 'Delimiter', '+'));
s = regexp(terms{1},'[(]?([0-9.]+)[)]? ([A-Za-z0-9_\-\[\]]+)','tokens');
for iterm = 1:length(s)
    if ~isempty(s{iterm})
        stoich = stmul*str2double(s{iterm}{1}{1});
        metab = s{iterm}{1}{2};
    else
        stoich = 1*stmul;
        metab = terms{1}{iterm};
    end                
    tf = strcmp(metab, model.metab);
    if any(tf)
        model.S(tf, irxn) = stoich;
    else
        model.S(imetab, irxn) = stoich;
        model.metab{imetab} = metab;
        imetab = imetab + 1;
    end
end
if ~isempty(poseff)
    for ipeff = 1:length(poseff)
        tfpos = strcmp(poseff{ipeff},model.metab);
        if any(tfpos)
           model.Z(irxn,tfpos) = 1;
        else
           model.Z(irxn,imetab) = 1;
           model.metab{imetab}= poseff{ipeff};
           imetab = imetab+1;
        end
    end
% else
%     model.Z(irxn,logical(model.S(:,irxn))) = 0;
end
if ~isempty(negeff)
    for ineff = 1:length(negeff)
        tfneg = strcmp(negeff{ineff},model.metab);
        if any(tfneg)
           model.Z(irxn,tfneg) = -1;
        else
           model.Z(irxn,imetab) = -1;
           model.metab{imetab}= negeff{ineff};
           imetab = imetab+1;
        end        
    end
% else
%     model.Z(irxn,logical(model.S(:,irxn))) = 0;
end
return

function [effector] = effectors(eff_string)
    %separate terms into a vector
    if ~isempty(eff_string)
        eff_string = strtrim(strrep(eff_string,'"',''));        
        commapos = strfind(eff_string,',');
        effector = cell(length(commapos)+1,1);
        if ~isempty(commapos)
            ipar = 1;
            while ipar <= length(commapos)
                if ipar == 1
                    effector{ipar} =...
                    strtrim(eff_string(1:commapos(ipar)-1));
                else
                    effector{ipar} =...
                    strtrim(eff_string(commapos(ipar-1)+1:commapos(ipar)-1));
                end
                ipar = ipar + 1;
            end
            effector{ipar} =...
            strtrim(eff_string(commapos(ipar-1)+1:end));
        else
            effector{1} = strtrim(eff_string(1:end));
        end
    else
        effector = [];          
    end
return
