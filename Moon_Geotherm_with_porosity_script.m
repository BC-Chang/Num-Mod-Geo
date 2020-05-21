H = 7.38E-12; % mean heat production [W/kg]
rho = 3300;   % mean density [kg/m^3]
R = 1738E3;   % Radius of the moon [m]
T0 = 250;     % mean surface temperature [K]
kappa = 3.3;  % mean thermal conductivity [W/(m K)]
phi0 = 0.25;  % surface porosity
hr = 30E3;    % decay depth
% Porosity correlation
phi = @(r) phi0*exp(-(R-r)/hr);

% Build Grid and operators
Grid.xmin = 0; 
Grid.xmax = R; 
Grid.Nx   = 300;
Grid.geom = 'spherical_r';
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
kappa_eff = (1-phi(Grid.xc))*kappa;
Kd = comp_mean(kappa_eff,-1,Grid);

L =-D*Kd*G; 
fs = (1-phi(Grid.xc))*rho*H;

% Build boundary conditions
BC.dof_dir   = Grid.dof_xmin;
BC.dof_f_dir = Grid.dof_f_xmin;
BC.g         = [T0];
BC.dof_neu   = []; 
BC.dof_f_neu = []; 
BC.qb        = [];
[B,N,fn] = build_bnd(BC,Grid,I);

% Solve for temperature and heat flux
u = solve_lbvp(L,fs+fn,B,BC.g,N);
q = comp_flux(D,Kd,G,u,fs,Grid,BC);

% Plot solution
subplot 411
plot(Grid.xc/1e3,fs*1e6,'b-','markerfacecolor','w','markersize',6)
ylabel '(1-\phi)\rho H [\muW/m3]'

subplot 412
plot(Grid.xc/1e3,kappa_eff,'b-','markerfacecolor','w','markersize',6)
ylabel '(1-\phi)\kappa [W/(m K)]'
ylim([2 4])

subplot 413
plot(Grid.xf/1e3,q*1e3,'b-','markerfacecolor','w','markersize',6), hold off
ylabel 'q [mW/m^2]'
ylim([0 15])

subplot 414
plot(Grid.xc/1e3,u,'b-','markerfacecolor','w','markersize',6)
xlabel 'r [km]', ylabel 'T [K]'
