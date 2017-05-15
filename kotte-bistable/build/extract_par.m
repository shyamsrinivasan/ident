function [parameter,lb,ub] = extract_par(par_string)
    %separate terms into a vector
    if ~isempty(par_string)
        par_string = strtrim(strrep(par_string,'"',''));        
        parm = strsplit(par_string,',');                      
        parameter = zeros(length(parm),1);
        lb = zeros(length(parm),1);
        ub = zeros(length(parm),1);
        if ~isempty(parm)
            ipar = 1;
            while ipar <= length(parm)
                parm{ipar} = strtrim(parm{ipar});
                bnd_op = strfind(parm{ipar},'[');
                bnd_cl = strfind(parm{ipar},']');
                if ~isempty(bnd_op) && ~isempty(bnd_cl)
                    %if there are bounds specified                    
                    Klub = strsplit(parm{ipar}(bnd_op+1:bnd_cl-1));                        
                    lb(ipar) = str2double(strtrim(Klub{1}));
                    ub(ipar) = str2double(strtrim(Klub{2}));
                else%No brackets
                    Kpar = str2double(parm{ipar});
                    if Kpar ~= 1
                        parameter(ipar) = Kpar;                            
                    end                            
                end
                ipar = ipar + 1;
            end
        end
    else
        parameter = [];
        lb = [];
        ub = [];
        %Default coefficients are defined within respective modules        
    end
end