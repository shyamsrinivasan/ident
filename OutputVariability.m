function sigma = OutputVariability(conc)
sigma = zeros(length(model.mets),length(model.mets));
for iout = 1:length(model.mets)
    for imet = 1:length(model.mets)
        outMean = mean(conc(iout,:));
        sigma(iout,imet) = sum((conc(iout,:)-outMean).^2)/(length(conc(iout,:))-1);
    end
end