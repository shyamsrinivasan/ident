function [defparval] = ExtractParameters(terms,defparval)

if nargin < 2
    defparval = struct();
end

nexpt = length(terms);
for exptnum = 1:nexpt
    compos = strfind(terms{exptnum},',');
    
    %defparval = struct();
    npar = 0;
    if ~isempty(compos)
        npar = length(compos);        
    end
    
    parameters = {'repcoeff';'accoeff';'rephill';'srate';'drate';'kmax';'pdrate'};
    
    npar = npar + 1;
    parterm = cell(npar,1);
    
    ipar = 1;
    while ipar <= npar
        
        if npar > 1 && ipar == npar
            parterm{ipar} = terms{exptnum}(compos(ipar-1)+1:end);
        elseif npar > 1 && ipar > 1
            parterm{ipar} = terms{exptnum}(compos(ipar-1)+1:compos(ipar)-1);
        elseif npar > 1 && ipar == 1
            parterm{ipar} = terms{exptnum}(1:compos(ipar)-1);
        else%npar == 1
            parterm{ipar} = terms{exptnum}(1:end);
        end
        
        eqpos = strfind(parterm{ipar},'=');
        parname = parterm{ipar}(1:eqpos-2);
        parval = parterm{ipar}(eqpos+2:end);
        
        if any(strcmpi(parname,parameters))
            org_val = defparval.(parameters{strcmpi(parname,parameters)});
            
            fprintf('Changing Default %s from %4.2g to %4.2g\n',...
                parameters{strcmpi(parname,parameters)},org_val,str2double(parval));
            
            defparval.(parameters{strcmpi(parname,parameters)}) =...
                                                        str2double(parval);
        end
        ipar = ipar + 1;
        
        
    end
end
end