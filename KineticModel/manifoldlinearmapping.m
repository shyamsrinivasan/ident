% curve_map_R2_R3.m, by Daniel Lyons (March 2014)
% projects points in x1-x2 plane to the plane spanned by the vectors W1,W2
% in the x1-x2-x3 coordinate systems
% Modifications include var and fun names and vectorization

function xnew = manifoldlinearmapping(x1,x2,W1,W2)

% perform linear mapping of unit circle onto plane in R3 spanned by W1,W2
% [W1 W2]; %transformation matrix

T = [W1 W2];
x1new = T(1,:)*[x1;x2];
x2new = T(2,:)*[x1;x2];
x3new = T(3,:)*[x1;x2];

xnew = [x1new;x2new;x3new];


% xnew = [W1 W2]*[x1;x2];


