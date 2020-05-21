%% Problem paramters
% Physical paramters
yr2s = 60^2*24*365.25; % [s/yr] convecrion from years to seconds
tmax = 15*1e4*yr2s;    % [s]    maximum simulation time
zmax = 5e3;            % [m]    Depth of simulation domain 
Gamma = 25/1000;       % [C/m]  geothermal gradient 
k = 1e-6;              % [m^2/s] thermal diffusivity 
v_e = -5e3/1e6/yr2s;   % erosion rate [m/s]

% Numerical paramters
Nt = 100;              % [-]    number of timesteps
dt = tmax/Nt;          % [s]    timestep size
theta = 0;             % [-]    timestepping paramter (0 = Backward Euler)

% Analytic expressions
T0 = @(z) Gamma*z;
Tb = @(t,v) Gamma*(zmax-v*t);  % Use the first part of analytic solution as bottom BC

%% Bulid grid and operators
% General discrete operators
Grid.xmin = 0; Grid.xmax = zmax; Grid.Nx = 100;
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
fs = sparse(Grid.Nx,1); % no heat generation

v = v_e*ones(Grid.Nfx,1);    % velocity vector
A = flux_upwind(v,Grid);    % upwind flux matrix

% Problem specific matrices
L = D*(A - k*G);      % Steady advection diffusion operator
Im = I + dt*(1-theta)*L;     % Implicit matrix
Ex = I - dt*theta*L;     % Explicit matrix

%% Boundary conditions
BC.dof_dir   = [Grid.dof_xmin; Grid.dof_xmax];
BC.dof_f_dir = [Grid.dof_f_xmin; Grid.dof_f_xmax];
BC.g = [0; 0]; % initialize to correct size, second entry has to be updated every timestep
BC.dof_neu = [];
BC.dof_f_neu = [];
BC.qb = [];
[B,N,fn] = build_bnd(BC,Grid,I);

%% Initial condition
u0 = T0(Grid.xc); % pre-established geotherm
u = u0;

%% timestepping loop
for i = 1:Nt
    time = dt*i;
    BC.g = [0; Tb(time,v_e)];
    u = solve_lbvp(Im,fs+Ex*u,B,BC.g,N);
    
    subplot 121
    plot(u0,Grid.xc/1e3,':',u,Grid.xc/1e3,'-')
    set(gca, 'YDir','reverse')
    xlabel('Temperature [C]')
    ylabel('Depth relative to current surface [km]')
    xlim([0 150]), ylim([0 5])
    
    subplot 122
    plot(u0,Grid.xc/1e3,':',u,(Grid.xc-v_e*time)/1e3)
    set(gca, 'YDir','reverse')
    xlabel('Temperature [C]')
    ylabel('Depth relative to initial surface [km]')
    xlim([0 150]), ylim([0 5])
    drawnow
end