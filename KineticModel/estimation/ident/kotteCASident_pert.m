function [ode,dfx_sym,D2FX,oderhs,x,p_all,ident_c,p_useless,acetate] =...
        kotteCASident_pert(idx,nc,nf,npert)

x = casadi.SX.sym('x',nc*npert,1);
p_all = casadi.SX.sym('p_all',12,1);
ident_c = casadi.SX.sym('ident_c',1,1);
p_useless = casadi.SX.sym('p_useless',3,1);
acetate = casadi.SX.sym('acetate',1,1);
% if idx==1
%     p = [ident_c;p_all];
% else
    p = [p_all(1:idx-1);ident_c;p_all(idx:end);p_useless;acetate];
% end
% p = [p_all(1:idx-1);ident_c;p_all(idx+1:end)];

flux = cell(6,1);

% parameters
K1ac = p(1);    % or 0.02
K3fdp = p(2);
L3fdp = p(3)*1e6;
K3pep = p(4);
K2pep = p(5);
vemax = p(6);        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = p(7);        % or 0.45
ne = p(8);             % or 2
d = p(9);
V4max = p(10);
k1cat = p(11);   
V3max = p(12);    
V2max = p(13);   
ac = p(17);

flux1 = k1cat.*x(3:nc:nc*npert).*ac./(ac+K1ac);
flux2 = vemax.*(1-1./(1+(KeFDP./x(2:nc:nc*npert)).^ne));
ratio = 1+x(2:nc:nc*npert)./K3fdp;
flux3 = V3max.*(ratio-1).*(ratio).^3./...
            (ratio.^4+L3fdp.*(1+x(1:nc:nc*npert)./K3pep).^(-4));
flux4 = V2max.*x(1:nc:nc*npert)./(x(1:nc:nc*npert)+K2pep);  
flux5 = V4max.*x(1:nc:nc*npert);
flux6 = d.*x(3:nc:nc*npert);

oderhs = [];
for i = 1:npert
    oderhs = [oderhs;flux1(i) - flux4(i) - flux5(i);...
                     flux4(i) - flux3(i);...
                     flux2(i) - flux6(i)];
end             
ode = casadi.Function('ode',{x,p_all,ident_c,p_useless,acetate},{oderhs}); 

dfx_sym = [];
D2FX = [];





