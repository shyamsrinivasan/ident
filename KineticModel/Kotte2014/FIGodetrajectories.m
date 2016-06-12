function [hfig,ha] =...
FIGodetrajectories(xdyn,ival,xeq,datatype,idx,hfig,ha)

end

function plot3Dtrajectories
xdata = xdyn(idx(1),:);
ydata = xdyn(idx(2),:);
zdata = xdyn(idx(3),:);

% plot trajectory
plot3(xdata,ydata,zdata)
hold on

% plot initial and final values
plot3(ival)
plot3(xeq)
end

function plot2Dtrajectories

end