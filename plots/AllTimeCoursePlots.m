function varargout = AllTimeCoursePlots(outsol,model,metid,flxid)
% wrapper around timecourseplots for plotting results from an ensemble of
% models
if nargin<4
    flxid = [];
end

nmodels = size(outsol,2);
hfig1 = [];
hfig2 = [];
ha1 = [];
ha2 = [];

% performance improvements for faster plots with lower memory burden
% combine all data from structure arrays into NaN separated matrices
npt = numel(outsol(1).t);
nmet = size(outsol(1).y,1);
nflx = size(outsol(1).flux,1);
allt = zeros(npt*nmodels+nmodels,1);
ally = zeros(nmet,npt*nmodels+nmodels);
allf = zeros(nflx,npt*nmodels+nmodels);
try
    allt(1:npt+1) = [outsol(1).t;NaN];
catch
    allt(1:npt+1) = [outsol(1).t NaN]';
end
ally(:,1:npt+1) = [outsol(1).y NaN(nmet,1)];
allf(:,1:npt+1) = [outsol(1).flux NaN(nflx,1)];
for im = 2:nmodels    
    try
        allt(((im-1)*(npt+1)+1):im*(npt+1)) = [outsol(im).t;NaN];
    catch
        allt(((im-1)*(npt+1)+1):im*(npt+1)) = [outsol(im).t NaN]';
    end
    ally(:,((im-1)*(npt+1)+1):im*(npt+1)) = [outsol(im).y NaN(nmet,1)];
    allf(:,((im-1)*(npt+1)+1):im*(npt+1)) = [outsol(im).flux NaN(nflx,1)];
end
if ~isempty(metid)
    [hfig1,ha1] = timecourseplots(allt,ally,...
                                  1,metid,model,hfig1,ha1);
end
if ~isempty(flxid)
    [hfig2,ha2] = timecourseplots(allt,allf,...
                                  2,flxid,model,hfig2,ha2);
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
