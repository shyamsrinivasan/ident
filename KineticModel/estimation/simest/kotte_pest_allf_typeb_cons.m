function cons =...
kotte_pest_allf_typeb_cons(xexp,vexp,nc,nf,npert,x,p,flux,vareps,p_usl,ac)

p_all = [p;p_usl;ac];

% parameters
K1ac = p_all(1);    % or 0.02
K3fdp = p_all(2);
L3fdp = p_all(3)*1e6;
K3pep = p_all(4);
K2pep = p_all(5);
vemax = p_all(6);        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = p_all(7);        % or 0.45
ne = p_all(8);             % or 2
d = p_all(9);
V4max = p_all(10);
k1cat = p_all(11);   
V3max = p_all(12);    
V2max = p_all(13);   
ac = p_all(17);

fx1 = k1cat.*x(3:nc:nc*npert).*ac - flux(1:nf:nf*npert).*(ac+K1ac);
fx2 = -vemax + (vemax - flux(2:nf:nf*npert)).*(1+(KeFDP./x(2:nc:nc*npert)).^ne);
ratio = 1+x(2:nc:nc*npert)./K3fdp;
fx3 = V3max.*(ratio-1).*(ratio).^3 -...
      flux(3:nf:nf*npert).*(ratio.^4+L3fdp.*(1+x(1:nc:nc*npert)./K3pep).^(-4));
fx4 = V2max.*x(1:nc:nc*npert) - flux(4:nf:nf*npert).*(x(1:nc:nc*npert)+K2pep);
fx5 = V4max.*x(1:nc:nc*npert) - flux(5:nf:nf*npert);
fx6 = d.*x(3:nc:nc*npert)-flux(6:nf:nf*npert);

% concentration noise cons
x_noise_lb = x - xexp*(1+vareps(1));
x_noise_ub = -x + xexp*(1-vareps(1));
      
% flux noise cons
v_noise_lb = flux - vexp*(1+vareps(2));
v_noise_ub = -flux + vexp*(1-vareps(2));      

cons = [];
for i = 1:npert
    cons = [cons;fx1(i);fx2(i);fx3(i);fx4(i);fx5(i);fx6(i)];
end
cons = [cons;...
        x_noise_lb;x_noise_ub;...
        v_noise_lb;v_noise_ub];





