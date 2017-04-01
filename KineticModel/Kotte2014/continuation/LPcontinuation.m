function LPcontinuation(x0,ap,funame,pvec)

global sys
sys.gui.pausespecial=1;  %Pause at special points 
sys.gui.pausenever=0;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

[x0,v0]=init_LP_LP(funame,x0,pvec',ap);
opt=contset;
opt=contset(opt,'MaxNumPoints',300);
opt=contset(opt,'MinStepSize',0.00001);
opt=contset(opt,'MaxStepSize',0.01);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Eigenvalues',1);
[x,v,s,h,f]=cont(@limitpoint,x0,v0,opt);

