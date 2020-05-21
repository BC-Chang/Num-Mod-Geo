%% Build grids and ops
mu = 1;
Gridp.xmin = -1; Gridp.xmax = 1; Gridp.Nx = 100;
Gridp.ymin =  0; Gridp.ymax = 1; Gridp.Ny =  50;
Grid = build_stokes_grid(Gridp);
[D,Edot,Dp,Gp,Z,I]=build_stokes_ops(Grid);
[Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc);

% Stokes operators
A = 2*mu*D*Edot; % 
L = [A, -Gp;...
     Dp, Z];
fs = spalloc(Grid.N,1,Grid.y.N);

%% Boundary conditions
BC.dof_dir = [Grid.dof_pene;...
              Grid.dof_pc];
BC.g = [zeros(Grid.N_pene,1);...
        0];
[B,N,fn] = build_bnd(BC,Grid,I);

%% Temperature field
T = zeros(Grid.p.N,1); T(Xc<0) = 1;
T = reshape(T,Grid.p.Ny,Grid.p.Nx);
% Averaged to the faces
T_f = comp_mean(T,1,1,Grid.p);
T_f = T_f(sub2ind(size(T_f),1:size(T_f,1),1:size(T_f,2)))';

fs(Grid.x.N+1:Grid.x.N+Grid.y.N) = T_f(Grid.x.N+1:end);

u = solve_lbvp(L,fs+fn,B,BC.g,N);
v = u(1:Grid.p.Nf); p = u(Grid.p.Nf+1:end);
PSI = comp_streamfun(v,Grid.p);

%% Plot
[Xp,Yp] = meshgrid(Grid.x.xc,Grid.y.yc);
psi_max = max(PSI(:));
colormap('bone')
contourf(Xc,Yc,reshape(T,Grid.p.Ny,Grid.p.Nx)), hold on
contour(Xp,Yp,PSI,psi_max*[.2:.2:1],'r-')
colorbar
axis equal
xlabel 'x', ylabel 'y'
legend('T','\Psi')
