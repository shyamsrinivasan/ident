function [pactnr,pactdr] = pregactivity(bindaff,act_ind,inh_ind,rule,data)
pactnr = 0;
pactdr = 0;
    if ~isempty(rule)
        orpos = strfind(rule,'|');
        andpos = strfind(rule,'&');
        if ~isempty(orpos) %If rule is made of ORs 
            if any(act_ind) && any(inh_ind)
                pactnr = sum(bindaff(act_ind));
                pactdr = sum(bindaff(act_ind))+sum(bindaff(inh_ind).^data.rephill);
            elseif any(act_ind)
                pactnr = sum(bindaff(act_ind));
                pactdr = sum(bindaff(act_ind));
            elseif any(inh_ind)
                pactnr = 1;
                pactdr = sum(bindaff(inh_ind).^data.rephill);
            end   
        end
        if ~isempty(andpos)%If rule is made of ANDs
            if any(act_ind) && any(inh_ind)
                pactnr = prod(bindaff(act_ind));
                pactdr = sum(bindaff(act_ind))+sum(bindaff(inh_ind).^data.rephill);
            elseif any(act_ind)
                pactnr = prod(bindaff(act_ind));
                pactdr = sum(bindaff(act_ind));
            elseif any(inh_ind)
                pactnr = 1;
                pactdr = sum(bindaff(inh_ind).^data.rephill);
            end           
        end
        %If it is a single regulator
        if isempty(orpos) && isempty(andpos)
            if any(act_ind)
                pactnr = sum(bindaff(act_ind));
                pactdr = sum(bindaff(act_ind));
            elseif any(inh_ind)
                pactnr = 1;
                pactdr = sum(bindaff(inh_ind).^data.rephill);
            end
        end
    end
return