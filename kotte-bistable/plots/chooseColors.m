function ColorSpec = chooseColors(ncolor,color)
if nargin < 2    
    yc = 0;
else
    yc = 1;
end
load('C:\Users\shyam\Documents\Courses\CHE1125Project\mat_files\Colors.mat');
if ~yc    
    Clrnd = floor(1 + (76-1)*rand(ncolor,1));
else
    if length(color)>1
        Clrnd = cellfun(@(x)strcmpi(Colors,x),color,'UniformOutput',false);
        Clrnd = cellfun(@(x)find(x),Clrnd,'UniformOutput',false);
        Clrnd = cell2mat(Clrnd);
    else
        Clrnd = repmat(find(strcmpi(Colors,color)),ncolor,1);
    end
end
ColorSpec = cell(length(Clrnd),1);
for j = 1:length(Clrnd)
    try 
        ColorSpec{j} = rgb(Colors{Clrnd(j)});
    catch
        ColorSpec{j} = rgb('Black');
    end
end

return