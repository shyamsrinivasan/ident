% runDynamicsims
% simulate from one or more limit points or within their boundaries 
% as initial conditions to test stability/dynamic behaviour of points
% load the relevant simulation dataset
load();

% needed variables: alliidpvec,alliidxeq,alliidfeq;
npts = size(alliidpvec,1);
% determine the #parameters/combinations that have been changed
ndp = size(alliidpvec,3);
for idp = 1:ndp
    % determine actual # parameters that have been changed
    npar = length(alliidpvec(1,:,idp));
    if npts>1
        diffpar = find(alliidpvec(1,:,idp)~=alliidovec(2,:,idp));
    end
end
