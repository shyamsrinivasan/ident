function obj = estimation_obj(x,pstruct)

if isfield(pstruct,'vnoisy')
    vnoisy = pstruct.vnoisy;
end

vmodel = convkin_kotte(x,pstruct);
obj = norm(vmodel-vnoisy);