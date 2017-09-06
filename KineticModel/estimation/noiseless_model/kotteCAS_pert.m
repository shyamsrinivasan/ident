function [ode,flux,oderhs,fluxeq,x,p_unch,dfx,dfxfun] = kotteCAS_pert(nc,npert)

x = casadi.SX.sym('x',nc*npert,1);
p_unch = casadi.SX.sym('p_unch',11,1);
p_pert = casadi.SX.sym('p_pert',3*npert,1);

% parameters
K1ac = p_unch(1);    % or 0.02
K3fdp = p_unch(2);
L3fdp = p_unch(3)*1e6;
K3pep = p_unch(4);
K2pep = p_unch(5);
vemax = p_unch(6);        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = p_unch(7);        % or 0.45
ne = p_unch(8);             % or 2
d = p_unch(9);
V4max = p_unch(10);
ac = p_unch(11);
k1cat = p_pert(1:3:3*npert);   
V3max = p_pert(2:3:3*npert);    
V2max = p_pert(3:3:3*npert);   

flux1 = k1cat(1:3).*x(3:nc:nc*npert).*ac./(ac+K1ac);
flux2 = vemax.*(1-1./(1+(KeFDP./x(2:nc:nc*npert)).^ne));
ratio = 1+x(2:nc:nc*npert)./K3fdp;
flux3 = V3max(1:3).*(ratio-1).*(ratio).^3./...
            (ratio.^4+L3fdp.*(1+x(1:nc:nc*npert)./K3pep).^(-4));
flux4 = V2max(1:3).*x(1:nc:nc*npert)./(x(1:nc:nc*npert)+K2pep);  
flux5 = V4max.*x(1:nc:nc*npert);
flux6 = d.*x(3:nc:nc*npert);

fluxeq = [];
for j = 1:npert
    fluxeq = [flux1(j);flux2(j);flux3(j);flux4(j);flux5(j);flux6(j)];
end

oderhs = [];
for i = 1:npert
    oderhs = [oderhs;flux1(i) - flux4(i) - flux5(i);...
                     flux4(i) - flux3(i);...
                     flux2(i) - flux6(i)];
end   
p = [p_unch;p_pert];

dfx = jacobian(oderhs,x);
dfxfun = casadi.Function('dfxfun',{x,p},{dfx});
ode = casadi.Function('ode',{x,p},{oderhs}); 
flux = casadi.Function('flux',{x,p},{fluxeq});                                        

                                        