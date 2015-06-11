%function [coeff,callflag] = coeffinit(regrules,defparval,callflag)
%**************************************************************************
%Initialize unknown activation and repression coefficients to default values

%Inputs
%regrules - List of regulators for which coefficients are to be assigned
%defvalueGene - Structure of default parameter values

%Optional Inputs
%callflag - Flag to indicate whether values will be written to file [{0}|1]

%Outputs
%coeff - A cell array of coefficeints and inhibition type

%October 2013 - version 1.0
%July 2014 - Version 2.0
%**************************************************************************
function [coeff,callflag] = coeffinit(regrules,defparval,callflag)
if nargout < 3
    callflag = 0;
end
if ~callflag
    %fprintf('Switching to Default Values \n');
    fprintf('Default Coefficients \n (1)Activation:%4.3d\n',defparval.accoeff);
    fprintf('(2)Repression:%4.3d\n (3)Dual:%4.3d\n',...
            defparval.repcoeff,0.00);        
    callflag = 1;
end
opbr = strfind(regrules,'(');
clbr = strfind(regrules,')');
nopbr = length(opbr);
orterms = {};
andterms = {};

coeff_br = cell(nopbr,1);
ntotal = 0;
ibr = 1;
while ibr <= nopbr
    orpos = strfind(regrules(opbr(ibr)+1:clbr(ibr)-1),'|');
    if ~isempty(orpos)
        orterms = strsplit(regrules(opbr(ibr)+1:clbr(ibr)-1),'|');
    end
    andpos = strfind(regrules(opbr(ibr)+1:clbr(ibr)-1),'&');
    if ~isempty(andpos)
        andterms = strsplit(regrules(opbr(ibr)+1:clbr(ibr)-1),'&');
    end
    if ~isemptyr(orterms) && ~isemptyr(andterms)
        coeff_rule = cell(length(orterms)+length(andterms),1);
    elseif ~isemptyr(orterms)
        coeff_rule = cell(length(orterms),1);
    elseif ~isemptyr(andterms)
        coeff_rule = cell(length(andterms),1);
    else
        orterms = strsplit(regrules(opbr(ibr)+1:clbr(ibr)-1),'|');
        coeff_rule = cell(length(orterms),1);
    end    
    kterm = 1;
    if ~isemptyr(orterms)
        jorterm = 1;        
        while jorterm <= length(orterms)
            [ormatch] = regexp(orterms{jorterm},'(\S+.?)(\W+.?)','tokens');             
            switch ormatch{1}{2}
                case '+'
                    coeff_rule{kterm}{1} = sprintf('%6.4g',defparval.accoeff);                   
                case '-'
                    coeff_rule{kterm}{1} = sprintf('%6.4g',defparval.repcoeff);
                case '+/-'
                    coeff_rule{kterm}{1} = sprintf('%6.4g',defparval.accoeff);
            end
            coeff_rule{kterm}{2} = ormatch{1}{2};
            kterm = kterm + 1;                        
            jorterm = jorterm + 1;
        end
    elseif ~isemptyr(andterms)
        jandterm = 1;
        while jandterm <= length(andterms)
            [andmatch] = regexp(andterms{jandterm},'(\w+.?)\[(\W+.?)\]','tokens');
            switch andmatch{1}{2}
                case '+'
                    coeff_rule{kterm}{1} = sprintf('%s',defparval.accoeff);                    
                case '-'
                    coeff_rule{kterm}{1} = sprintf('%s',defparval.repcoeff);              
                case '+/-'
                    coeff_rule{kterm}{1} = sprintf('%s',defparval.accoeff);
            end
            coeff_rule{kterm}{2} = andmatch{1}{2}; 
            kterm = kterm + 1;                                   
            jandterm = jandterm + 1;
        end
    end 
    coeff_rule = coeff_rule(~cellfun('isempty',coeff_rule));
    coeff_br{ibr} = coeff_rule;
    ntotal = ntotal + length(coeff_br{ibr});
    ibr = ibr + 1;    
end
coeff = cell(ntotal,1);
ktot = 1;
for kbr = 1:nopbr
    for jterm = 1:length(coeff_br{kbr})
        coeff{ktot} = coeff_br{kbr}{jterm};
        ktot = ktot + 1;
    end
end
    
end
