function [hfig,ha] = plotrndivals(rndivals,ssid,alleq,id,datatype,hfig,ha)
if nargin<5
    datatype = 2;
end
if nargin<6
    hfig = [];
end
if nargin<7
    ha = [];
end

colorSpec = chooseColors(2,{'Green','Purple'});

% collect all steady states
ss1ival = rndivals(ssid==1,:);
ss1eq = alleq(:,ssid==1);
Point1.MarkerFaceColor = colorSpec{2};
Point1.MarkerEdgeColor = colorSpec{2};
ss2ival = rndivals(ssid==2,:);
ss2eq = alleq(:,ssid==2);
Point2.MarkerFaceColor = colorSpec{1};
Point2.MarkerEdgeColor = colorSpec{1};

if ~isempty(ss1ival) && ~isempty(ss1eq)
    [hfig,ha] =...
    FIGmssEqIvalPerturbations(ss1ival',ss1eq(:,1),datatype,id,hfig,ha,Point1);  
    drawnow
end

% for iss1 = 1:size(ss1ival,1)
%     [hfig,ha] =...
%     FIGmssEqIvalPerturbations(ss1ival(iss1,:)',ss1eq(:,iss1),datatype,id,hfig,ha,Point1);  
%     drawnow
% end
if ~isempty(ss2ival) && ~isempty(ss2eq)
    [hfig,ha] =...
    FIGmssEqIvalPerturbations(ss2ival',ss2eq(:,1),datatype,id,hfig,ha,Point2);  
    drawnow
end

% for iss2 = 1:size(ss2ival,1)
%     
% end
