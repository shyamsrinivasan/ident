load('C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\Kinetic Model\N2m.mat');
bounds.vl = zeros(FBAmodel.nt_rxn,1);
% bounds.vl(bounds.vl==0) = -1;
bounds.vl(5) = -1;
bounds.vu = zeros(FBAmodel.nt_rxn,1);
bounds.vu(bounds.vu==0) = 1;
bounds.Vuptake = 1;
[vLPmax,vLPmin] = solveLP(FBAmodel,'P','P5',bounds,2)