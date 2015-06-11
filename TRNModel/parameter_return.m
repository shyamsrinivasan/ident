%function [coefficients,srates] = return_parameters(temp_param)
%return individual parameters matrices from a combined vector of parameters
function [final_coeff,brate] = parameter_return(temp_par,model)
    ngene = length(model.Gene);
    tnreg = length(model.Regulators);
    
    [gene,prot] = find(model.RS);
    coeff_mat = sparse(gene,prot,1,ngene,tnreg);
    %coeff_mat = sparse(logical(model.RS),logical(model.RS),1,ngene,tnreg);
    brate = zeros(ngene,1);

    currpos = 0;
    for igene = 1:ngene         
        nreg = length(find(model.RS(igene,:)));
        brate(igene,1) = temp_par(currpos+1);
        coeff = temp_par(currpos+2:currpos+1+nreg);     
        currpos = currpos+1+nreg;
        ireg = 1;
        while ireg <= nreg        
            tf = find(model.RS(igene,:));
            coeff_mat(igene,tf(ireg)) = coeff(ireg);            
            ireg = ireg+1;
        end
    end
    final_coeff = sparse(coeff_mat);    
end