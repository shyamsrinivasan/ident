function varargout = AllTimeCoursePlots(outsol,model,metid,flxid)
if nargin<4
    flxid = [];
end
% wrapper around timecourseplots for plotting results from an ensemble of
% models

nmodels = size(outsol,2);
hfig1 = [];
hfig2 = [];
ha1 = [];
ha2 = [];
for im = 1:nmodels
    if ~isempty(metid)
        [hfig1,ha1] = timecourseplots(outsol(im).t,outsol(im).y,...
                                      1,metid,model,hfig1,ha1);
    end
    if ~isempty(flxid)
        [hfig2,ha2] = timecourseplots(outsol(im).t,outsol(im).y,...
                                      2,flxid,model,hfig2,ha2);
    end
end

if nargout==1
    varargout{1} = hfig1;
elseif nargout==2
    varargout{1} = hfig1;
    varargout{2} = ha1;
elseif nargout==3
    varargout{1} = hfig1;
    varargout{2} = ha1;
    varargout{3} = hfig2;
elseif nargout==4
    varargout{1} = hfig1;
    varargout{2} = ha1;
    varargout{3} = hfig2;
    varargout{4} = ha2;
end
