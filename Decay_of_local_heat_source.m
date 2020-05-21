%% Physical properties
kappa = 4;
cp = 2e3;
rho = 1e3;
Tini = 300;
T0 = 150;

s2yr = 60^2*24*365;
tmax = 2e6*s2yr;
Nt = 5e1;
dt = tmax/Nt;
theta = 0;

%% Build grid and ops
Grid.xmin = 0;
Grid.xmax = 100e3;
Grid.Nx = 5e2;
Grid.geom = 'cylindrical_r';
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
L = -D*kappa*G; M = I*rho*cp;

IM = @(theta,dt) M + (1-theta)*dt*L;
EX = @(theta,dt) M - theta*dt*L;
fs = 0;
Im = IM(theta,dt); Ex = EX(theta,dt);

%% Boundary condtions
Param.dof_dir = Grid.dof_xmax;
Param.dof_f_dir = Grid.dof_f_xmax;
Param.g = [T0];
Param.dof_neu = [];
Param.dof_f_neu = [];
Param.qb = [];
[B,N,fn] = build_bnd(Param,Grid,I);

%% Initial condition
u0 = ones(Grid.N,1)*T0;
r1 = 4e3; r2 = 12e3;
a = find(Grid.xc >= r1 & Grid.xc <=r2);
u0(1:a(1)) = Tini;
u0(a) = Tini - (Tini-T0)/(r2-r1)*(Grid.xc(a)-r1);


plot(Grid.xc/1e3,u0), hold on
ylim([100 300])

%% Time evolution
time = 0; u = u0;
for i = 1:Nt
    u = solve_lbvp(Im,Ex*u,B,Param.g,I);
    plot(Grid.xc/1e3,u,'r-')
    ylim([100 300])
    drawnow
end