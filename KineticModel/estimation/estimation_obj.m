function obj = estimation_obj(x,pstruct)
% vmodel - flux from convinience kientics model
% vnoisy - flux from 

if isfield(pstruct,'vnoisy')
    vnoisy = pstruct.vnoisy;
end
if isfield(pstruct,'flxh')
    flxh = pstruct.flxh;
end
if isfield(pstruct,'c')
    c = pstruct.c;
end

vmodel = flxh(c,x);
obj = norm(vmodel-vnoisy);