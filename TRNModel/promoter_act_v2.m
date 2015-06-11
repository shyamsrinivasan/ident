function [pactnr,pactdr,psrate_prot,pact] =...
         promoter_act_v2(terms,srate,bindaff,defparval,RS,Protein)
%**************************************************************************
%Function to calculate single promoter activity in conjunction with
%singlepromoteractivity.m
%******Inputs(s)
%terms
%srate
%bindaff
%defparval
%RS
%Protein
%******Outputs(s)
%pactnr
%pactdr
%psrate_prot
%pact
%
%April 01 2014
% - Created to work similar to v5 of singlepromoteractivity but for all
%   cases of boolean logic. Called in v7 of singlepromoteractivity

%**************************************************************************


nprots = length(terms);
old_tf = zeros(length(Protein),nprots);
iprot = 1;
tfcount = 1;
psrate_prot = zeros(nprots,1);
use_srate = zeros(nprots,1);
logic = {};

while iprot <= nprots
    tf = strcmp(terms{iprot}{1},Protein);
    %=========First Protein/TF in the Gene Regulatory Rule=========
    if iprot == 1             
        switch RS(1,tf)
            case 1
                pactnr = bindaff(tfcount);                
                pactdr = bindaff(tfcount);
                use_srate(iprot) = 1;
            case -1                         
                pactnr = 1;
                pactdr = (bindaff(tfcount))^defparval.rephill;  
                                
            case 2 %pactivity of dual regulators is same as +ve TFs
                pactnr = bindaff(tfcount);%                         
                pactdr = bindaff(tfcount);
                use_srate(iprot) = 1;
        end  
        psrate_prot(iprot) = srate(1,tf);
        old_tf(:,iprot) = tf;       
        tfcount = tfcount + 1;
        if length(terms{iprot}) > 1 && isempty(logic)
            logic = terms{iprot}{2};        
        end                
    elseif iprot > 1 && iprot < nprots 
        %==============All but first or last regulator=============
        [pactnr,pactdr,psrate_prot] =...
         promoter_act(logic,RS,iprot,pactnr,pactdr,psrate_prot);
        old_tf(:,iprot) = tf; 
        if isempty(logic)
            logic = terms{iprot}{2};   
        end
    else %iprot == nprots
        %========================Last Regulator====================
        [pactnr,pactdr,psrate_prot] =...
         promoter_act(logic,RS,iprot,pactnr,pactdr,psrate_prot);
        old_tf(:,iprot) = tf;        
    end 
    iprot = iprot + 1;

end

if pactnr == 1
    pact = sum(psrate_prot)*pactnr/(length(psrate_prot)*(1+pactdr));
else
    pact = pactnr/(1+pactdr);
end



function [pactnr,pactdr,psrate] =...
         promoter_act(logic,RS,kprot,pactnr,pactdr,psrate)
    %prev_effect = RS(1,logical(old_tf(:,kprot-1)));
    if logic == '|'
        switch RS(1,tf)
            case 1
                if pactnr == 1
                    pactnr = pactnr*bindaff(tfcount);
                else
                    pactnr = pactnr + bindaff(tfcount);
                end
                pactdr = pactdr + bindaff(tfcount);
                psrate(kprot) = srate(1,tf);            
                
            case -1
                pactnr = pactnr*1;
                pactdr = pactdr + (bindaff(tfcount))^defparval.rephill; 
                psrate(kprot) = srate(1,tf);
                
            case 2
                if pactnr == 1
                    pactnr = pactnr*bindaff(tfcount);
                else
                    pactnr = pactnr + bindaff(tfcount);
                end
                pactdr = pactdr + bindaff(tfcount); 
                psrate(kprot) = srate(1,tf);                        
        end       
    elseif logic == '&'
        switch RS(1,tf)
            case 1
                if pactnr == 1
                    pactnr = pactnr*bindaff(tfcount);
                    use_srate(kprot) = 1;
                else
                    if any(use_srate == 1)
                        
                        pactnr = pactnr*bindaff(tfcount);
                        
                        use_srate(use_srate==1) = 2;
                        use_srate(kprot) = 1;
                        
                    elseif any(use_srate == 2)
                        
                        pactnr = pactnr*bindaff(tfcount);
                        
                    end
                    %test how many srates have been added?
                    %How to test?
                    %Step assign srates to psrates
                    %If psrate is used as psrate(kprot-1) change
                    %psrate(krpot-1) to 0 => psrate(krpot-1) has been
                    %used or use a flag variable
                    %flag = 0 -  Not used
                    %flag = 1 - Used once
                    %flag = 2 - Used in an average?
                end        
                pactdr = pactdr + bindaff(tfcount);
                psrate(kprot) = srate(1,tf);
                
            case -1                
                pactnr = pactnr*1;
                pactdr =...
                pactdr + pactdr*(bindaff(tfcount))^defparval.rephill;
            
                psrate(kprot) = srate(1,tf);
                
            case 2
                %test how many srates have been added?
                    %How to test?
                    %Step assign srates to psrates
                    %If psrate is used as psrate(kprot-1) change
                    %psrate(krpot-1) to 0 => psrate(krpot-1) is has been
                    %used or use a flag variable
                    %flag = 0 -  Not used
                    %flag = 1 - Used once
                    %flag = 2 - Used in an average?
                if pactnr == 1
                    pactnr = pactnr*bindaff(tfcount);
                    use_srate(kprot) = 1;
                else
                    if any(use_srate == 1)
                        
                        pactnr = pactnr*bindaff(tfcount);
                    
                        use_srate(use_rate==1) = 2;
                        use_srate(kprot) = 1;
                        
                    elseif any(use_srate == 2)
                        
                        pactnr = pactnr*bindaff(tfcount);
                        
                    end                    
                end        
                pactdr = pactdr + bindaff(tfcount);
                psrate(krpot) = srate(1,tf);                           
        end        
    end 
    tfcount = tfcount + 1;
end

end