%% Set problem paramters
h = 40e3;      % crustal thickness [m]
hr = 9e3;     % decay depth [m]
rho = 3000;    % crustal density [kg/m^3]
H0 = 1e-9;     % surface heat production [W/m^3]
kappa = 3.5;  % thermal conductivity of crust [W/m/K]
T0 = 15;     % mean surface temperature [C]
qm = 18e-3;     % Mantle heat flow [W/m^2]

%% Build grid and operators
Grid.xmin = 0; 
Grid.xmax = h; 
Grid.Nx = 100;
Grid.geom = 'cartesian';
Grid = build_grid(Grid);
[D,G,I]=build_ops(Grid);
L = -D*kappa*G;
fs = rho*H0*exp(-Grid.xc/hr);      % r.h.s. - radiogenic heat production

%% Set boundary conditions
BC.dof_dir   = Grid.dof_xmin;         % identify cells on Dirichlet bnd
BC.dof_f_dir = Grid.dof_f_xmin;         % identify faces on Dirichlet bnd
BC.dof_neu   = Grid.dof_xmax;         % identify cells on Neumann bnd
BC.dof_f_neu = Grid.dof_f_xmax;         % identify faces on Neumann bnd
BC.g  = [T0];                % set bnd temperature
BC.qb = [qm];                % set bnd heat flux
[B,N,fn] = build_bnd(BC,Grid,I);  % Build constraint matrix and basis for its nullspace

%% Solve linear boundary value problem & plot solution
u = solve_lbvp(L,fs+fn,B,BC.g,N);
q = comp_flux(D,kappa,G,u,fs,Grid,BC);

subplot 121
plot(u,Grid.xc/1e3)
xlabel 'T [K]', ylabel 'z [km]'
subplot 122
plot(q*1e3,Grid.xf/1e3)
xlabel 'q [W/m^2]', ylabel 'z [km]'
