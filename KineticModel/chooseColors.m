function ColorSpec = chooseColors(ncolor)
load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Colors.mat');
Clrnd = floor(1 + (76-1)*rand(ncolor,1));
ColorSpec = cell(length(Clrnd),1);
for j = 1:length(Clrnd)
    try 
        ColorSpec{j} = rgb(Colors{Clrnd(j)});
    catch
        ColorSpec{j} = rgb('Black');
    end
end
return