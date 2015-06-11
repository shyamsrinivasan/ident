%calculate eigen vectors and eigen values from covariance matrices
covariance = cov(conc');
[eigenvector,eigenvalue] = eig(covariance);
