% curve_map_R2_R3.m, by Daniel Lyons (March 2014)
% projects points in x1-x2 plane to the plane spanned by the vectors W1,W2
% in the x1-x2-x3 coordinate systems
% Modifications include var and fun names and vectorization

function xnew = manifoldlinearmapping(x1,x2,W1,W2)

% perform linear mapping of unit circle onto plane in R3 spanned by W1,W2
% [W1 W2]; %transformation matrix

% x1_new = zeros(1,npoints+1);
% x2_new = zeros(1,npoints+1);
% x3_new = zeros(1,npoints+1);

xnew = [W1 W2]*[x1;x2];


% for ii=1:npoints+1
%     x_mat = [x1(ii); x2(ii)];
%     x1_new(ii) = T(1,:)*x_mat;
%     x2_new(ii) = T(2,:)*x_mat;
%     x3_new(ii) = T(3,:)*x_mat;
% end


% [x1new,x2_new,x3_new] =...
% manifoldprojection(x1,x2,eigvec_stable(:,1),eigvec_stable(:,2),npoints);