% hobj = get(ha1,'Children');
set(0,'CurrentFigure',hf1);
set(hf1,'CurrentAxes',ha1);
iter = 39;
pid = 1;
id = 79;
oldid = 79;
while id >= 0    
    xdata = get(hobj(id),'XData');
    ydata = get(hobj(id),'YData');    
    if oldid-id == 4    
        pid = pid+1;
        oldid = id;
    end
    ptid = ['P' num2str(pid)];
    text(xdata,ydata,ptid);
    iter = iter-1;    
    id = 2*iter+1;
end