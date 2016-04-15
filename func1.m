function F = func1(z)
%rhs matrix
func = @(x)(x.^2-2*x+1);
x = [0;0.155051;0.644949;1];
K = 3;
L = zeros(K+1,1);
Ldot = zeros(K+1,1);
a = zeros(K+1,1);
z0 = -3;
for jk = 1:K
    for j = 0:K        
        [L(j+1),Ldot(j+1)] = Lagrange(x(jk+1),j,x);
        a(j+1) = Ldot(j+1);
    end
    A(jk,:) = a';    
end

F = A(:,1)*z0 + A(:,2:end)*z'-func(z');