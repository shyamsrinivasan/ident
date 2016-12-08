function Manifold2DPlot(x,y,z)

% delaunay triangulation of surface
t = Delaunay2_5D([x' y' z']);

% get indices of triangles inside boundary box and outside
cut_val = 6; 
[tri_inside,tri_boundary,tri_outside] = tri_calc(x,y,z,t,cut_val);

% redraw boundary triangles so they end at the boundary value "cut_val"
[tri_inside_new,tri_boundary_saved,x,y,z] = polydraw(x,y,z,tri_boundary,tri_inside,cut_val);

figure
hold on
h6=trisurf(tri_inside_new,x,y,z,'facecolor','interp',...
                                'edgecolor','none',...
                                'edgelighting','phong',...
                                'facelighting','phong');
                            
set(h6,'facecolor',[0 1 0]);%green
set(h6,'facealpha',0.25);  
xlabel('pep a.u.','fontsize',18);ylabel('fdp a.u.','fontsize',18);
zlabel('E a.u.','fontsize',18);
set(gca,'fontsize',18);
