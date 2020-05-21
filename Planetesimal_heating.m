yr2s = 60^2*24*365.25; % [s/yr] number of seconds per year

% Primary parameters (from Hevey and Sanders 2006, Table 1)
Param.k   = 2.1;            % [W/(m K)] Thermal conductivity 
Param.cp  = 837;            % [J/(kg K)] specific heat capacity 
Param.rho = 3400;           % [kg/m3] density 
Param.lam = 9.5e-7/yr2s;    % [1/s] decay constant converted from 1/yr to 1/s
Param.H0  = 1.9e-7;         % [W/kg] initial power per unit mass
Param.T0  = 250;            % [K] initial and surface temperature
Param.R   = 100e3;          % [m] Planetesimal radius 
Param.t0  = 2.25e6*yr2s;    % [s] accretion time

% Derived parameters
Param.kappa = Param.k/(Param.rho*Param.cp);                        % [m/s2] thermal diffusivity 
Param.A0    = @(t0) Param.rho*Param.H0*exp(-Param.lam*t0);   % [W/m3] initial power per unit volume at accretion

% Analytic solution
%tau_cen = linspace(0,150e6,1e3)'*yr2s;
%T_cen = PlanetesimalAnalytic(Param,tau_cen,1,1e3);

Nt = 500;
tmax = 150e6*yr2s;
dt = tmax/Nt;

% Build grid
Grid.xmin = 0;
Grid.xmax = Param.R;
Grid.Nx = 1e2;
Grid.geom = 'spherical_r';
Grid = build_grid(Grid);

% Discrete Operators
[D,G,I] = build_ops(Grid);
L = -D*Param.k*G;
M = I*Param.rho*Param.cp;
IM = @(theta,dt) M + (1-theta)*dt*L;
EX = @(theta,dt) M - theta*dt*L;
a = ones(Grid.Nx,1);
fs_int = @(tnew,told) -a.*Param.A0(Param.t0)/Param.lam*(exp(-Param.lam*tnew) - exp(-Param.lam*told))/dt;
% Build BC's
BC.dof_dir = Grid.dof_xmax;
BC.dof_f_dir = Grid.dof_f_xmax;
BC.g = [Param.T0];
BC.dof_neu = [];
BC.dof_f_neu = [];
BC.qb = [];
[B,N,fn] = build_bnd(BC,Grid,I);

% Time stepping loop
theta = 0;

u = Param.T0*ones(Grid.Nx,1); ucen = zeros(Nt,1); time = ucen;
for i=1:Nt
    time(i) = dt*i;
    told = (i-1)*dt; tnew = time(i);
    u = solve_lbvp(IM(theta,dt),EX(theta,dt)*u+dt*(fs_int(tnew,told)+fn),B,BC.g,N);
    ucen(i) = u(1); 
end
figure
%plot((tau_cen+Param.t0)/yr2s/1e6,T_cen,'-'), hold on
plot((time+Param.t0)/yr2s/1e6,ucen,'--')
xlabel 'time [Ma]', ylabel 'T [K]'
xlim([0 150])
