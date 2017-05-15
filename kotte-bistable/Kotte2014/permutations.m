% permutations
p = nchoosek([0.01 0.1 1 10 100],3);
v = zeros(3,0);
for ip = 1:size(p,1)
    vp = perms(p(ip,:));
    v = [v vp'];
end
v = unique(v','rows');