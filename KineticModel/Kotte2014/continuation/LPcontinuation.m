function data = LPcontinuation(x0,ap,funame,pvec)

global sys
sys.gui.pausespecial=1;  %Pause at special points 
sys.gui.pausenever=0;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

% limit point conitnuation (2 free variables in ap)
[x0,v0]=init_LP_LP(funame,x0,pvec',ap);
opt=contset;
opt=contset(opt,'MaxNumPoints',300);
opt=contset(opt,'MinStepSize',0.00001);
opt=contset(opt,'MaxStepSize',0.01);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Eigenvalues',1);
[x,v,s,h,f]=cont(@limitpoint,x0,v0,opt);

% separation of variable and parameter vectors
if ~isempty(s1)    
    y = x1(1:length(x0),:);
    p = x1(length(x0)+1:end,:); % 2 free vairables
else
    y = [];
    p = [];
end

if ~isempty(s)
    data.s = s;
    data.x = x;
    data.f = f;
    data.v = v;
    data.h = h;
%      if ~isempty(flux)
%          data.flux = flux;
%      else
%          data.flux = [];
%      end
else
    data = [];
end
