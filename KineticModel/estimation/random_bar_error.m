a=[0,1,0,0;
4,3,2,1;
2,2,1,3;
1,0,0,0];
b=[0,1,0,0;
1,2,1,1;
1,1,1,2;
1,0,0,0];
ctrs = 1:4;
data = a;
figure(1)
hBar = bar(ctrs, data);
for k1 = 1:size(a,1)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(k1,:) = hBar(k1).YData;
end
hold on
errorbar(ctr, ydt, b, '.r')
hold off