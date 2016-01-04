function v = columnVector(v)


if size(v,1)>1 && size(v,2)==1
    v = v';
end