function obj = obj_flux1_noisy(x,id,val)

obj = sparse(id,1,val,7,1)'*x;