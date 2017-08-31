function [obj,var,par,x,p,flux,vareps,p_usl,ac,wts,ss_obj] =...
        kotte_pest_allf_typeb_obj(xexp,vexp,nc,nf,npert)
x = casadi.SX.sym('x',nc*npert,1);
p = casadi.SX.sym('p',13,1);
p_usl = casadi.SX.sym('p_usl',3,1);
ac = casadi.SX.sym('ac',1,1);
flux = casadi.SX.sym('flux',nf*npert,1);
vareps = casadi.SX.sym('vareps',2,1);
wts = casadi.SX.sym('wts',4,1);

p_all = [p;p_usl;ac];

% parameters
% K1ac = p_all(1);    % or 0.02
% K3fdp = p_all(2);
% L3fdp = p_all(3)*1e6;
% K3pep = p_all(4);
% K2pep = p_all(5);
% vemax = p_all(6);        % for bifurcation analysis: 0.7:0.1:1.3
% KeFDP = p_all(7);        % or 0.45
% ne = p_all(8);             % or 2
% d = p_all(9);
% V4max = p_all(10);
% k1cat = p_all(11);   
% V3max = p_all(12);    
% V2max = p_all(13);   
ac = p_all(17);

vmodel = flux;
v_error = vexp-vmodel;
v_norm = .5*dot(v_error,v_error);
x_error = xexp-x;
x_norm = .5*dot(x_error,x_error);

ss_obj = [];
for i = 1:npert
    flux_pert = flux(nf*(i-1)+1:nf*i);
    nlerhs = [flux_pert(1)-flux_pert(4)-flux_pert(5);...
              flux_pert(4)-flux_pert(3);...
              flux_pert(2)-flux_pert(6)];  
    ss_obj = [ss_obj;nlerhs];
end

ss_obj_norm = .5*dot(ss_obj,ss_obj);

var = [x;p;flux;vareps];
par = [p_usl;ac;wts];

obj = 10000*ss_obj_norm + wts(1)*v_norm + wts(2)*x_norm + wts(3)*vareps(1) + wts(4)*vareps(2);


