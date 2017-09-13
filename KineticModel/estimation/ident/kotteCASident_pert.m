function [ode,flux,D2FX,oderhs,x,p_var,p_ident,p_fixed,acetate] =...
        kotteCASident_pert(idx,nc,nf,npert)

x = casadi.SX.sym('x',nc*npert,1);
p_var = casadi.SX.sym('p_var',8,1);
p_ident = casadi.SX.sym('p_ident',1,1);
p_fixed = casadi.SX.sym('p_fixed',4,1);
acetate = casadi.SX.sym('acetate',1,1);
% rearrange parameters
p = [p_var(1:idx-1);p_ident;p_var(idx:end);p_fixed;acetate];

% parameters
K1ac = p(1);    % or 0.02
K3fdp = p(2);
L3fdp = p(3)*1e6;
K3pep = p(4);
K2pep = p(5);
V4max = p(6);
k1cat = p(7);   
V3max = p(8);    
V2max = p(9);
vemax = p(10);        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = p(11);        % or 0.45
ne = p(12);             % or 2
d = p(13); 
ac = p(14);

flux1 = k1cat.*x(3:nc:nc*npert).*ac./(ac+K1ac);
flux2 = vemax.*(1-1./(1+(KeFDP./x(2:nc:nc*npert)).^ne));
ratio = 1+x(2:nc:nc*npert)./K3fdp;
flux3 = V3max.*(ratio-1).*(ratio).^3./...
            (ratio.^4+L3fdp.*(1+x(1:nc:nc*npert)./K3pep).^(-4));
flux4 = V2max.*x(1:nc:nc*npert)./(x(1:nc:nc*npert)+K2pep);  
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
ode = casadi.Function('ode',{x,p_var,p_ident,p_fixed,acetate},{oderhs}); 
flux = casadi.Function('flux',{x,p_var,p_ident,p_fixed,acetate},...
                                {fluxeq});

% dfx_sym = [];
D2FX = [];





