clear
clc

%% Build Stokes Grid Function Call
mu = 1; % viscosity - value does not matter here

% Definition of pressure grid
Gridp.xmin = 0; Gridp.xmax = 1; Gridp.Nx = 5;
Gridp.ymin = 0; Gridp.ymax = 1; Gridp.Ny = 5;

Grid = build_stokes_grid(Gridp);

%% Build Stokes Operator Function Call
mu = 1; % viscosity - value does not matter here

% Definition of pressure grid
Gridp.xmin = 0; Gridp.xmax = 1; Gridp.Nx = 5;
Gridp.ymin = 0; Gridp.ymax = 1; Gridp.Ny = 5;

Grid = build_stokes_grid(Gridp);
[D,Edot,Dp,Gp,Z,I] = build_stokes_ops(Grid);

A = 2*mu*D*Edot; % 
L = [A, -Gp;...
     Dp, Z];

subplot 141
spy(D), title 'D'
subplot 142
spy(Edot), title 'Edot'
subplot 143
spy(A), title 'A'
subplot 144
spy(L), title 'L'