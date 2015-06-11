%Log transform data & add noise
function [logdata,varargout] = NoiseTransform(Data,Options)

%Log transforming data

Data(Data~=0) = log(Data(Data~=0));
logdata = Data;

%Generate Noise
NoiseSD = Options.NoiseSD;
ntimepoints = Options.MaxDataPoints;
nvars = size(Data,1);
Noise = randn(nvars,ntimepoints).*NoiseSD;

%Add Noise;
LogNoisyData = logdata + Noise;
%npoints = solverP.MaxDataPoints;
%nvars = size(Solution.initSS{2},1);
% mean = zeros(nvars,1);
% mean(mean==0) = 1e-9;
% Noise = randn(nvars,npoints).*5e-1;

varargout{1} = LogNoisyData;
varargout{2} = Noise;


end