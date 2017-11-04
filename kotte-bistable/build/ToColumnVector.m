function vec = ToColumnVector(vec)
[nr,nc] =size(vec);
if nc==1 && nr>nc
    vec = vec';
end