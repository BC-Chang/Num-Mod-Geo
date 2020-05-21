mu = 1;

%% Build staggered grids
Gridp.xmin = 0; Gridp.xmax = 1; Gridp.Nx = 50;
Gridp.ymin = 0; Gridp.ymax = 1; Gridp.Ny = 50;
Grid = build_stokes_grid(Gridp);

%% Build Stokes operators
[D,Edot,Dp,Gp,Z,I] = build_stokes_ops(Grid);
A = 2*mu*D*Edot; 
L = [A, -Gp; ...
    Dp, Z];
fs = spalloc(Grid.N,1,0);

%% BC's
BC.dof_dir = [Grid.dof_ymin_vx(2:end-1);... % tangential velocity on the top (ymin)
              Grid.dof_pene;...    % no penetration on all bnd's
              Grid.dof_pc];        % Pressure constraint
          
BC.g       = [ones(Grid.p.Nx-1,1);...       % tangential velocity on the top
              zeros(Grid.N_pene,1);...   % no penetration on all bnd's
              0];                           % Pressure constraint
[B,N,fn] = build_bnd(BC,Grid,I);

%% Solve for Stokes flow
u = solve_lbvp(L,fs+fn,B,BC.g,N);
v = u(1:Grid.p.Nf); 
p = u(Grid.p.Nf+1:end);
PSI = comp_streamfun(v,Grid.p);

%% Plot solution
[Xp,Yp] = meshgrid(Grid.x.xc,Grid.y.yc);

contour(Xp,Yp,PSI,10,'k')
set(gca,'ydir','reverse')
axis square
