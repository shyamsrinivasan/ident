% run reginVolume dynamic simulations from sampled code
load('regionVolumeSamplevals.mat');
% integrate
options = [];
[ppival,npival,allxeq,ssid,allfeq] =...
getPequilibrium(rndivals(20001:30000,:),model,pvec,options,opts,tspanf);
