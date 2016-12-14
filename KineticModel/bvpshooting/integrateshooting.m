function [tf,yf,status] = integrateshooting(fh,ti,tf,yi,opts)
status = 1;

[tdyn,ydyn] = ode45(fh,ti:0.1:tf,yi,opts);
tfnew = tdyn(end);
yf = ydyn(end,:)';

if tfnew ~= tf
    status = -1;
    tf = tfnew;
end