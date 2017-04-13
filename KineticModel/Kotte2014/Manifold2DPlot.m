function hfig = Manifold2DPlot(x,y,z,hfig,type,icolor)
if nargin<5
    type = 1;
end
if nargin<6
    icolor = [0 1 0]; % green
end
if nargin<4 || isempty(hfig)
    hfig = figure;
end
set(0,'CurrentFigure',hfig);
ha = gca;
set(ha,'NextPlot','add');

if type == 1
    % delaunay triangulation of surface
    t = Delaunay2_5D([x' y' z']);

    % get indices of triangles inside boundary box and outside
    cut_val = 2.5; 
    [tri_inside,tri_boundary,tri_outside] = tri_calc(x,y,z,t,cut_val);

    % redraw boundary triangles so they end at the boundary value "cut_val"
    [tri_inside_new,tri_boundary_saved,x,y,z] =...
    polydraw(x,y,z,tri_boundary,tri_inside,cut_val);

    hsurf=trisurf(tri_inside_new,x,y,z,'facecolor','interp',...
                                    'edgecolor','none',...
                                    'edgelighting','phong',...
                                    'facelighting','phong');
    % set(hsurf,'EdgeColor',icolor);
    set(hsurf,'FaceColor',icolor);
    set(hsurf,'FaceAlpha',0.25); 
elseif type == 2
    plot3(x,y,z,'LineStyle','none','Marker','.','Color',[0 .5 0]);
end
[xlabel,ylabel,zlabel] = getKotteaxislabels(3,2,[1 2 3]);
setKotteproperties(3,ha,xlabel,ylabel,zlabel);

% xlabel('pep a.u.','fontsize',18);ylabel('fdp a.u.','fontsize',18);
% zlabel('E a.u.','fontsize',18);
% set(gca,'fontsize',18);
