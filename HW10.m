%% This example is just a variation on the heat flow around an inclusion so we have nice heat-flow streamlines
Kin = 100;   % Thermal conductivity in interior
rad_max = 2; % Major half axis of ellipse
rad_min = .2; % Minor half axis of ellipse
theta = 45;   % rotation of the ellipse

% Numerical grid
Lx = 10; Ly = 5;
Grid.xmin = -Lx/2; Grid.xmax = Lx/2;  Grid.Nx = 4*100;
Grid.ymin = -Ly/2; Grid.ymax = Ly/2;  Grid.Ny =  4*50;

Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
fs = sparse(Grid.N,1);
[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
%% Thermal conductivity field
K = ones(Grid.Ny,Grid.Nx);
% boundary of ellipse
t = linspace(0,2*pi,1000);
x_ellipse = rad_max*cos(t)*cos(deg2rad(theta)) - rad_min*sin(t)*sin(deg2rad(theta));
y_ellipse = rad_max*cos(t)*sin(deg2rad(theta)) + rad_min*sin(t)*cos(deg2rad(theta));

in = inpolygon(Xc(:),Yc(:),x_ellipse,y_ellipse);
dof_in  = Grid.dof(in==1);
dof_out = Grid.dof(in==0);

K(dof_in) = Kin;
Kd = comp_mean(K,-1,1,Grid);
L = -D*Kd*G; 

% Boundary conditions
BC.dof_dir = [Grid.dof_xmin;   Grid.dof_xmax]; 
BC.dof_f_dir = [Grid.dof_f_xmin;   Grid.dof_f_xmax]; 
BC.dof_neu = []; BC.dof_f_neu = [];
BC.g = [ones(Grid.Ny,1);...
     zeros(Grid.Ny,1)];
BC.qb = [];

% Build boundary operators
[B,N,fn] = build_bnd(BC,Grid,I);

% Solve system
h = solve_lbvp(L,fs+fn,B,BC.g,N); 
q = comp_flux(D,Kd,G,h,fs,Grid,BC);

%% Compute streamfunction - t
[Xp,Yp] = meshgrid(Grid.xf,Grid.yf);
[PSI, psi_min1, psi_max1] = comp_streamfun(q,Grid);


%% Plotting
Ns = 2*11; Nh = 2*Ns;
[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])

colormap('cool')
subplot 211
contourf(Xc,Yc,K), colorbar, hold on
xlabel('x'), ylabel('y')
axis equal
set(gca,'ytick',[-5:5])
title('Thermal conductivity')

subplot 212
hmin = min(h); hmax = max(h);
psi_min = min(PSI(:)); psi_max = max(PSI(:));

h_cont = linspace(hmin,hmax,Nh);
psi_cont = linspace(psi_min,psi_max,Ns);

[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
contour(Xc,Yc,reshape(h,Grid.Ny,Grid.Nx),h_cont,'b-','linewidth',2), hold on
[Xp,Yp] = meshgrid(Grid.xf,Grid.yf);
contour(Xp,Yp,PSI,psi_cont,'r','linewidth',2)
xlim([Grid.xmin Grid.xmax]), ylim([Grid.ymin Grid.ymax])
xlabel 'x [m]', ylabel 'y [m]'
plot(x_ellipse,y_ellipse,'k-');
axis equal
set(gca,'ytick',[-5:5])
xlabel 'x [m]', ylabel 'y [m]'
title 'Temperature and Heatflow'
legend('T','q')