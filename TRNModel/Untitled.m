%function to perform sensitivity analysis to parameters
%function [] = sensAnalysis(parName,parScale)
parList = {'Kcat','K','KIact','KIihb','ptrate','Kb','Kub','Coefficient','brate','gmax'};
%           1      2   3       4       5        6    7     8             9
%           10
Ns = size(parName,1);%Rows of parameters
YS = zeros(nvar,Ns);

nsampl = 10;
%Find and Sample Desired Parameters
for ip = 1:Ns
    %Change one parameter at a time and simulate model
    %Change Parameter    
    tf_par = strcmpi(parName{ip,1},parList);
    if any(tf_par)%Parameter Exists
        pVal = pScale(ip,1) + (pScale(ip,2)-pScale(ip,1))*betarnd(1.5,4.5,1,nsampl);       
        %Simulate Model for revised parameters 
        for jsampl = 1:nsampl
            %Simulate model for each sample of the parameter value
            [model.Coefficient,model.brate] = parameter_return(data.par,model);
            pmeter = model.pmeter;            
            tf_gen = strcmpi(parName{ip,2},model.Gene);
            tf_prot = strcmpi(parName{ip,3},model.Enzyme);
            tf_reg = strcmpi(parName{ip,3},model.Regulators);
            tf_met = strcmpi(parName{ip,4},model.Metabolites);             
            %Select the switched parameter
            if find(tf_par) == 1
                model.Kcat(:) = pVal(jsampl);
            elseif find(tf_par) >=2 && find(tf_par) <= 4              
                if any(tf_prot) && and(tf_met)                    
                    pmeter.(parName{ip,1})(tf_met,tf_prot) = pVal(jsampl);
                end               
            elseif find(tf_par) > 4
                switch parName{ip,1}                             
                    case 'ptrate'
                        if any(tf_reg) && any(tf_gen)
                            model.ptrate(tf_reg,tf_gen) = pVal(jsampl);
                        end
                    case 'Kb'
                        if any(tf_reg)
                            model.Kub(tf_reg) = pVal(jsampl);
                        end
                    case 'Kub'
                        if any(tf_reg)
                            model.Kub(tf_reg) = pVal(jsampl);
                        end
                    case 'Coefficient'
                        if any(tf_gen) && any(tf_reg)                        
                            model.(parName{ip,1})(tf_gen,tf_reg) = pVal(jsampl);
                        end                    
                    case 'brate'
                        if any(tf_gen)
                            model.brate(tf_gen) = pVal(jsampl);
                        end
                    case 'gmax'
                        model.gmax = pVal(jsampl);                    
                end                   
            end
            model.allpar = parameter_vector(model,length(model.Gene));
            model.pmeter = pmeter;
        end        
    else
        fprintf('Parameter(s) %s is non-Existent',parName{ip,1});
    end
end
return




