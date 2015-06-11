%function to sample parameters
function pVal = MCsample_parameters(parName,pScale,nsampl)
parList = {'Kcat','K','KIact','KIihb','ptrate','Kb','Kub','Coefficient','brate','gmax'};
%           1      2   3       4       5        6    7     8             9
%           10
Ns = size(parName,1);%Rows of parameters
pVal = zeros(Ns,nsampl);
%Find and Sample Desired Parameters
if ~isempty(parName)
    for ip = 1:Ns
        %Change one parameter at a time and simulate model
        %Change Parameter    
        tf_par = strcmpi(parName{ip,1},parList);
        if any(tf_par)%Parameter Exists
            pVal(ip,:) = pScale(ip,1) + (pScale(ip,2)-pScale(ip,1))*betarnd(1.5,4.5,1,nsampl);             
        else
            fprintf('Parameter(s) %s is non-Existent',parName{ip,1});
        end
    end
end            
return




