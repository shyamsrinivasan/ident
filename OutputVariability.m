function sigma = OutputVariability(val)
nvar = size(val,1);
sigma = zeros(nvar,1);
for iout = 1:nvar
    outMean = mean(val(iout,:));
    sigma(iout) = sum((val(iout,:)-outMean).^2)/(length(val(iout,:))-1);
%     for imet = 1:length(model.mets)
%         
%     end
end