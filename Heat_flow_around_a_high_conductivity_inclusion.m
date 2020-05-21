Kin = 100;   % Thermal conductivity in interior
rad_max = 1; % Major half axis of ellipse
rad_min = 1; % Minor half axis of ellipse
theta = 0;   % rotation of the ellipse

% Numerical grid
Lx = 10; Ly = 5;
Grid.xmin = -Lx/2; Grid.xmax = Lx/2;  Grid.Nx = 4*100;
Grid.ymin = -Ly/2; Grid.ymax = Ly/2;  Grid.Ny =  4*50;

Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
fs = sparse(Grid.N,1);

%% Thermal conductivity field
K = ones(Grid.Ny,Grid.Nx); % initilize to unity background
% boundary of ellipse
t = linspace(0,2*pi,1000);
x_ellipse = rad_max*cos(t)*cos(deg2rad(theta)) - rad_min*sin(t)*sin(deg2rad(theta));
y_ellipse = rad_max*cos(t)*sin(deg2rad(theta)) + rad_min*sin(t)*cos(deg2rad(theta));

[xc,yc] = meshgrid(Grid.xc,Grid.yc);
in = inpolygon(xc,yc,x_ellipse,y_ellipse);
dof_in = find(in);  % dof's inside the ellipse
dof_out = find(~in); % dof's outside ellipse

K(dof_in) = Kin; %update K inside inclusion
Kd     = comp_mean(K,-1,1,Grid);  % conductivity matrix
L      = -D*Kd*G;  % steady head conduction operator

% Boundary conditions
BC.dof_dir   = [Grid.dof_xmin; Grid.dof_xmax];
BC.dof_f_dir = [Grid.dof_f_xmin; Grid.dof_f_xmax];
BC.g = [ones(length(Grid.dof_xmin),1); zeros(length(Grid.dof_xmax),1)];
BC.dof_neu   = [];
BC.dof_f_neu = [];
BC.qb        = [];
[B,N,fn]     = build_bnd(BC,Grid,I);
% Solve system
u = solve_lbvp(L,fs,B,BC.g,N); % solve for temperature field
q = comp_flux(D,Kd,G,u,fs,Grid,BC); % compute heat flow vector

%% Compute streamfunction (we will discuss this in two weeks - helps to plot heat flow)
[Xp,Yp] = meshgrid(Grid.xf,Grid.yf);
[PSI, psi_min1, psi_max1] = comp_streamfun(q,Grid);


%% Plotting
Ns = 2*11; Nh = 2*Ns;
[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
colormap('cool')
subplot 211
contourf(Xc,Yc,K), colorbar, hold on
xlabel('x'), ylabel('y')
axis equal
set(gca,'ytick',[-5:5])
title('Thermal conductivity')

subplot 212
plot_flownet_class(Nh,Ns,u,PSI,'b-','r',Grid), hold on
plot(x_ellipse,y_ellipse,'k-');
axis equal
set(gca,'ytick',[-5:5])
xlabel 'x [m]', ylabel 'y [m]'
title 'Temperature and Heatflow'
legend('T','q')
