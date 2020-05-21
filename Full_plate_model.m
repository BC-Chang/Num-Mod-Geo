%% Physical parameters
yr2s = 60^2*24*365.25; % [s/yr]
kappa = 3.0;          % [W/(m K)] thermal coductivity
cp = 1.171e3;         % [J/(kg K)] specific heat capacity 
rho = 3300;           % [kg/m^3] density
k = kappa/(rho*cp);   % [m^2/s] thermal diffusivity
Ts = 0;               % [C] surface temperature
Tr = 1450;            % [C] basal temperature
tmax = 200e6*yr2s;    % [s] maximum age of the seafloor
zmax = 95e3;          % [m] Plate thickness
v_plate  = 2e-2/yr2s; % [m/s] plate velocity [m/s]
xmax = v_plate*tmax;

%% Build grid and operators
Grid.xmin = 0; Grid.xmax = xmax; Grid.Nx = 75;
Grid.ymin = 0; Grid.ymax = zmax; Grid.Ny = 75;
Grid = build_grid(Grid); 
[D,G,I] = build_ops(Grid);
fs = sparse(Grid.Nx*Grid.Ny,1);

% Velocity field
vx = v_plate*ones(Grid.Nfx,1);
vz = zeros(Grid.Nfy,1);
v = [vx; vz];
A = flux_upwind(v,Grid);
L = D*(A-k*G); % steady advection diffusion operator

%% Build BC's
% Dirichlet
BC.dof_dir   = [Grid.dof_ymin(2:end); Grid.dof_xmin(2:end-1); Grid.dof_ymax; Grid.dof_xmin(1)];
BC.dof_f_dir = [Grid.dof_f_ymin(2:end); Grid.dof_f_xmin(2:end-1); Grid.dof_f_ymax; Grid.dof_f_xmin(1)];
BC.g         = [Ts*ones(Grid.Nx-1,1); Tr*ones(Grid.Ny-2,1); Tr*ones(Grid.Nx,1); (Ts+Tr)/2];

% Sort Dirichlet BC's
[BC.dof_dir,index] = sort(BC.dof_dir);
BC.dof_f_dir = BC.dof_f_dir(index);
BC.g = BC.g(index);

% Neumanm
BC.dof_neu   = [];
BC.dof_f_neu = [];
BC.qb        = [];
[B,N,fn]     = build_bnd(BC,Grid,I);

%% Solve for temperature field
u = solve_lbvp(L,fs,B,BC.g,N);
u(1:150)


%% Plotting
[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
Age_Myrs = Xc/v_plate/1e6/yr2s;
subplot 121
contourf(Xc/1e3,Yc/1e3,reshape(u,Grid.Ny,Grid.Nx),[100:100:1400]), hold on
colorbar('location','northoutside')
set(gca,'ydir','reverse','linewidth',2)
ylabel 'Depth [km]'
xlabel 'x [km]'
pbaspect([1 .4 1])
xlim([0 xmax]/1e3), ylim([0 zmax]/1e3)

subplot 122
contourf(Age_Myrs,Yc/1e3,reshape(u,Grid.Ny,Grid.Nx),[100:100:1400]), hold on
colorbar('location','northoutside')
set(gca,'ydir','reverse','linewidth',2)
ylabel 'Depth [km]'
xlabel 'age [Myr]'
pbaspect([1 .4 1])
xlim([0 tmax]/yr2s/1e6),ylim([0 zmax]/1e3)
