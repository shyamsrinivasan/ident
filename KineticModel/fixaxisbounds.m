function [minY,maxY] = fixaxisbounds(minY,maxY)
if maxY > 0
    if minY >= 1e-4 && maxY >= 1e-4
        if maxY-minY <= 1e-2
            maxY = maxY+1e-1;
            if minY-1e-1 <= 0
                minY = 0;
            else
                minY = minY-1e-1;
            end
        end
    else
        if maxY-minY <= 1e-4
            maxY = maxY+1e-4;
            if minY-1e-4 <= 0
                minY = 0;
            else
                minY = minY-1e-4;
            end
        end
    end
else
    tmin = -minY;
    tmax = -maxY;
    if tmin >= 1e-4 && tmax >= 1e-4
        if tmin-tmax <= 1e-2
            tmin = tmin+1e-1;
            if tmax-1e-1 <= 0
                tmax = 0;
            else
                tmax = tmax-1e-1;
            end
        end
    else
        if tmin-tmax <= 1e-4
            tmin = tmin+1e-4;
            if tmax-1e-4 <= 0
                tmax = 0;
            else
                tmax = tmax-1e-4;
            end
        end
    end
    maxY = -tmax;
    minY = -tmin;
end
return