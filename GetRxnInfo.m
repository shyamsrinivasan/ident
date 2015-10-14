function [model,reverse,imetab,irxn] = GetRxnInfo(model,imetab,irxn,parstr,rxnstring)

ipt = 0;
%Building S and K matrices - separate terms into a vector
[par,Klb,Kub] = extract_par(parstr);    
    
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

%compartmentalization
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

%reaction RHS
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

if size(model.mets,2)>size(model.mets,1)
    model.mets = model.mets';
end
% model.nt_metab = size(model.S,1);
