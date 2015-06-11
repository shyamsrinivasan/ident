fgrow = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1];
freg = [20 16 12 10 8 6 4 2 1 .9 .6 .3 .1 .01 0];
rst = zeros(length(fgrow),length(freg));
for igrow = 1:length(fgrow)
    for ireg = 1:length(freg)
        rst(igrow,ireg) = fgrow(igrow)*(1/50 + freg(ireg));
    end
end

for ireg = 1:length(freg)
    plot(fgrow,rst(:,ireg));
    hold on
end















