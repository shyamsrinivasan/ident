function xnew = manifoldlinearmapping(x1,x2,W1,W2)
% Adapted from Lyons (2014) curve_map_R2_R3.m
% project 2-D points onto 3-D plane spanned by vectors W1 and W2 

% perform linear mapping of unit circle onto plane in R3 spanned by W1,W2
% [W1 W2]; %transformation matrix

T = [W1 W2];
x1new = T(1,:)*[x1;x2];
x2new = T(2,:)*[x1;x2];
x3new = T(3,:)*[x1;x2];

xnew = [x1new;x2new;x3new];


% xnew = [W1 W2]*[x1;x2];


