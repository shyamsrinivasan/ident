function rr = algebraicSystem(Y,data,model,FBAmodel)
%Algebraic system at steady state
mdecay = 0.0031;
pdecay = 3.83e-6;
gr_flux = 0.8/3600;
alphac = 0.06;
alpha = 40/(233/(gr_flux*3600)^2+78);
beta = 60/(82.5/(gr_flux*3600)+148);
RNAP = 3e-5;
b = 5000;
k = 1e-4;%PAR
c1 = ((0.06*beta/k)/(pdecay+gr_flux));
c2 = (mdecay + gr_flux)/(RNAP);
c3 = alpha/(alphac*b);
rr = c1*c2*Y^2+(c2 - alpha/b*c1 - alphac*c1)*Y - alpha/b;
end