function [] = noEnzFlux()
VnoEnz = model.VnoEnz;
f_noEnz = zeros(length(VnoEnz,1));
for irxn = 1:length(VnoEnz)
    subsind = model.S(:,VnoEnz(irxn)) < 0;%substrate
    prodind = model.S(:,VnoEnz(irxn)) > 0;
    f_noEnz(irxn) = k*conecntration of substrate;
end
%Separate binding events from otehr spontaneous reactions
return