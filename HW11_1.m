yr2s = 60^2*24*365.25; % [s/yr]
mu = 1;
kappa = 3.0;        % [W/(m K)] thermal coductivity
cp = 1.171e3;       % [J/(kg K)] specific heat capacity 
rho = 3300;         % [kg/m^3] density
k = kappa/(rho*cp);  % [m^2/s] thermal diffusivity
Ts = 0;
Tr = 1200;
vs = 4e-1/yr2s;

% characteristic parameters
x_c = k/vs;

%% Build staggered grids for flow calculation
W = 200;
Gridp.xmin = 0; Gridp.xmax = W; Gridp.Nx = 1*W;
Gridp.ymin = 0; Gridp.ymax = W; Gridp.Ny = 1*W;
Grid_flow = build_stokes_grid(Gridp);
[Xx,Yx] = meshgrid(Grid_flow.x.xc,Grid_flow.x.yc);  % Middle of x-faces (location of vx on flow grid)
[Xy,Yy] = meshgrid(Grid_flow.y.xc,Grid_flow.y.yc);  % Middle of y-faces (location of vy on flow grid)
%% Build Stokes operators for flow calculation
[D_flow,Edot,Dp,Gp,Z,I_flow] = build_stokes_ops(Grid_flow);
A_flow = 2*mu*D_flow*Edot; % 
L_flow = [A_flow, -Gp;...
          Dp, Z];
fs_flow = spalloc(Grid_flow.N,1,0);

%% Build Flow BC's
BC.flow.dof_dir = [Grid_flow.dof_ymin_vt(2:end-1);... % tangential velocity on the 'top' (ymin) bnd 
              Grid_flow.dof_xmax_vt(2:end-1);...      % tangential velocity on the right (xmax) bnd
              Grid_flow.dof_ymax_vt(2:end-1);...      % tangential velocity on the 'bot' (ymax) bnd
              Grid_flow.dof_xmin_vt(2:end-1);...      % tangential velocity on the left (xmin) bnd
              Grid_flow.dof_pene;...                  % no penetration on all bnd's
              Grid_flow.dof_pc];                      % Pressure constraint
          
BC.flow.g  = [ones(Grid_flow.x.Nx-4,1);...       % tangential velocity on the 'top' (ymin) bnd
              ones(Grid_flow.y.Ny-4,1);...       % tangential velocity on the right (xmax) bnd
              -ones(Grid_flow.x.Nx-4,1);...       % tangential velocity on the 'bot' (ymax) bnd
              -ones(Grid_flow.y.Ny-4,1);...       % tangential velocity on the left (xmin) bnd
              zeros(Grid_flow.N_pene,1);...      % no penetration on all bnd's
              0];                                % Pressure constraint
          
[B_flow,N_flow,fn_flow] = build_bnd(BC.flow,Grid_flow,I_flow);

%% Solve for Stokes flow
u = solve_lbvp(L_flow,fs_flow+fn_flow,B_flow,BC.flow.g,N_flow);
v = u(1:Grid_flow.p.Nf); 
p = u(Grid_flow.p.Nf+1:end);

%% Plot flow solution
PSI = comp_streamfun(v,Grid_flow.p);
[Xp,Yp] = meshgrid(Grid_flow.x.xc,Grid_flow.y.yc);

contour(Xp,Yp,PSI,15,'k'), hold on
plot([0 W/2 W/2 0 0],[0 0 W/2 W/2 0],'r-','linewidth',1.5)
xlabel 'x', ylabel 'y'
axis equal tight

%% Transfer the velocity field to the transport grid
vx = v(1:Grid_flow.x.N);   % vx vector 
VX = reshape(vx,Grid_flow.x.Ny,Grid_flow.x.Nx);   % vx reshaped as matrix matching [Xx,Yx]
vy = v(Grid_flow.x.N+1:end);   % vy vector 
VY = reshape(vy,Grid_flow.y.Ny,Grid_flow.y.Nx);   % vy reshaped as matrix matching [Xy,Yy]

VX_trans = VX(1:round(Grid_flow.x.Ny/2),1:round(Grid_flow.x.Nx/2)); % VX matrix restricted to transport grid
VY_trans = VY(1:round(Grid_flow.y.Ny/2),1:round(Grid_flow.y.Nx/2)); % VY matrix restricted to transport grid
v_trans = [VX_trans(:); VY_trans(:)];  % v vector on transport grid

%% Solve Temperature field
% Transport grid
Grid_trans.xmin = 0; Grid_trans.xmax = W/2; Grid_trans.Nx = Gridp.Nx/2;
Grid_trans.ymin = 0; Grid_trans.ymax = W/2; Grid_trans.Ny = Gridp.Ny/2;
[Grid_trans] = build_grid(Grid_trans);
[Xc,Yc] = meshgrid(Grid_trans.xc,Grid_trans.yc); % location of cell centers in transport grid

[D_trans,G_trans,I_trans] = build_ops(Grid_trans);
A = flux_upwind(v_trans,Grid_trans);
L_trans = D_trans*(A-G_trans);
fs_trans = spalloc(Grid_trans.N,1,0);

%% Build Transport BC's
BC_trans.dof_dir   = [Grid_trans.dof_ymin; Grid_trans.dof_ymax];
BC_trans.dof_f_dir = [Grid_trans.dof_f_ymin; Grid_trans.dof_f_ymax];
BC_trans.g =         [ zeros(Grid_trans.Nx,1); ones(Grid_trans.Nx,1)];
 
BC_trans.dof_neu = [];
BC_trans.dof_f_neu = [];
BC_trans.qb = [];
[B_trans,N_trans,fn_trans] = build_bnd(BC_trans,Grid_trans,I_trans);

% Solve for advection-diffusion
T = solve_lbvp(L_trans,fs_trans+fn_trans,B_trans,BC_trans.g,N_trans);

%% Plot solution
figure;
subplot 121
contourf(Xc,Yc,reshape(T,Grid_trans.Ny,Grid_trans.Nx),50,'LineColor','none'), hold on
colorbar
contour(Xp,Yp,PSI,15,'k'), hold on
xlabel 'x_D', ylabel 'y_D'
title 'Dimensionless Temperature'
caxis([0 1])
xlim([0 100])
ylim([0 100])
set(gca,'ydir','reverse')
pbaspect([2 1 1])

subplot 122
Tdim = Ts + (Tr-Ts)*T;
contourf(Xc*x_c/1e3,Yc*x_c/1e3,reshape(Tdim,Grid_trans.Ny,Grid_trans.Nx),50,'LineColor','none'), hold on
colorbar
contour(Xp*x_c/1e3,Yp*x_c/1e3,PSI,15,'k'), hold on
caxis([Ts Tr])
xlabel 'x [km]', ylabel 'depth [km]'
title 'Temperature [C]'
xlim([0 100]*x_c/1e3)
ylim([0 100]*x_c/1e3)
set(gca,'ydir','reverse')
pbaspect([2 1 1])
