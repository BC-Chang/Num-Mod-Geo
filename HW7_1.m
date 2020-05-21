%% Grid and discrete operators
Grid.xmin = 0;  Grid.xmax = 2; Grid.Nx = 10;
Grid.ymin = -1; Grid.ymax = 1; Grid.Ny = 15;
Grid = build_grid(Grid); [Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
[D,G,I] = build_ops(Grid);
L = -D*G; fs = zeros(Grid.N,1);
flux = @(h) -G*h;
res = @(h,cell) L(cell,:)*h - fs(cell); 


%% Boundary conditions
Param.dof_dir   = [1; Grid.dof_xmin(2:Grid.Ny); Grid.dof_ymin(2:Grid.Nx)];
Param.dof_f_dir = [Grid.dof_f_xmin;Grid.dof_f_ymin(2:Grid.Nx)];
Param.g         = [.5;ones(Grid.Ny-1,1);zeros(Grid.Nx-1,1)];
Param.dof_neu   = [];
Param.dof_f_neu = [];
Param.qb        = [];
[B,N,fn] = build_bnd(Param,Grid,I);

%% Solve problem and compute fluxes
u = solve_lbvp(L,fs+fn,B,Param.g,N);
v = comp_flux_gen(flux,res,u,Grid,Param);


%% Compute advection operator
A = flux_upwind(v,Grid);

%% Plot results
[X2Y,Y2X] = interp_mat(Grid); % Matrices that interpolate fluxes from x to y faces and vice versa
[Xx,Yx] = meshgrid(Grid.xf,Grid.yc);
Vx = reshape(v(1:Grid.Nfx),Grid.Ny,Grid.Nx+1); 
Vy = reshape(Y2X*v(Grid.Nfx+1:Grid.Nf),Grid.Ny,Grid.Nx+1); 

contour(Xc,Yc,reshape(u,Grid.Ny,Grid.Nx),10), hold on
axis equal 
xlim([Grid.xmin Grid.xmax])
ylim([Grid.ymin Grid.ymax])
xlabel 'x', ylabel 'y'

quiver(Xx,Yx,1.5*Vx,1.5*Vy,'k')