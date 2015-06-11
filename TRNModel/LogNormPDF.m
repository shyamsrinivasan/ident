%function for determination of log-likelihood
%Liklehood
%L(Y|Theta) = PI(1/sqrt(2pi)sigma)*EXP(-(Data(i,j) - Model(i,j))/2sigma)
%Variables for likelihood for gaussian variable
%Expt value - y(t)
%Model value or mean - x(t)
%Variance - sigma
%NoisyData = NoisyData + randn(NumOfSpecies, D).*NoiseSD;

%Output = sum( -ones(D,1)*(0.5*log(2*pi*Variance)) - ((Values-Means).^2)./(2*(ones(D,1)*Variance))
%Solution.initSS{1} = timepoints
%Solution.initSS{2} = datamatrix

function [Output] = LogNormPDF( Values, Means, Variance)
%Output = Log Likelihood Estimate 


% Values:     D by 1 Vector
% Means:      D by 1 Vector
% Variance:   D by 1 Matrix

if size(Values, 2) > 1
    Values = Values';
end

if size(Means, 2) > 1
    Means = Means';
end

D = length(Values);


Output = sum( -ones(D,1)*(0.5*log(2*pi*Variance)) - ((Values-Means).^2)./(2*(ones(D,1)*Variance)) );

end


