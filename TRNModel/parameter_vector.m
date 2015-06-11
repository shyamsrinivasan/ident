%function temp_param = parameter_vector(coefficient,srate)
%Function to convert coefficients & srates into a single vector of
%parameters
%Order of parameters :[brate;Coefficients];
%1.brate    
%2.Coefficients    
function temp_par = parameter_vector(model,ngene)
    temp_par = [];      
    for jgene = 1:ngene
        coeff = model.Coefficient(jgene,:);
        coeff = coeff(coeff>0)';
        temp_par = [temp_par;model.brate(jgene);coeff];             
    end
end