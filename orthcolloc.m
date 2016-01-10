%Algorithm for solution to tiff ODEs using orthgonal collaction on finite
%elements

%N,i - number of elements
%K,j or k - number of collocation points within element
%pol - polynomial name or type whose roots form the points
%h - size of each element

%inputs
%x0 - intial point of integration
%xf - final point of integration
%x - total integration horizon

%z

N = 1;
K = 3;
pol = 'Gauss-Radau';

x = xf-x0;
h = x/N;

for iel = 1:N
    for jcp = 1:K
        
    end
end
