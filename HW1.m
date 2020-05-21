clear all
close all
clc

%% Problem 1
Grid.xmin = 0;
Grid.xmax = 1;
Grid.Nx = 10;
[Grid] = build_grid(Grid);

%% Problem 2
Grid.xmin = 0; Grid.xmax = 3; Grid.Nx = 30;
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid); L = D*G;
xa = linspace(Grid.xmin,Grid.xmax,3e2);
f = @(x) exp(cos(2*pi*x));
dfdx = @(x) -exp(cos(2*pi*x)).*sin(2*pi*x)*2*pi; 
d2fdx2 = @(x) -2*pi^2*exp(cos(2*pi*x)).*(2*cos(2*pi*x)+cos(4*pi*x)-1);

subplot 311
plot(xa,f(xa),'r',Grid.xc,f(Grid.xc),'bo')
xlabel 'x'
ylabel 'f'
legend('analytical','numerical')

subplot 312
plot(xa,dfdx(xa),'r',Grid.xf,G*f(Grid.xc),'bo')
xlabel 'x'
ylabel 'df/dx'
legend('analytical','numerical')

subplot 313
plot(xa,d2fdx2(xa),'r',Grid.xc,L*f(Grid.xc),'bo')
xlabel 'x'
ylabel 'd^2f/dx^2'
legend('analytical','numerical')

function [Grid] = build_grid(Grid)
% function [Grid] = build_grid(Grid)
% Author: Marc Hesse
% Date: 09/12/2014
% author: Bernard Chang
% date: 01/31/2020
% Description:
% This function computes takes in minimal definition of the computational
% domain and grid and computes all containing all pertinent information 
% about the grid. 
% Input:
% Grid.xmin = left boundary of the domain
% Grid.xmax = right bondary of the domain
% Grid.Nx   = number of grid cells
% 
% Output: (suggestions)
% Grid.Lx = scalar length of the domain
% Grid.dx = scalar cell width
% Grid.Nfx = number of fluxes in x-direction
% Grid.xc = Nx by 1 column vector of cell center locations
% Grid.xf = Nfx by 1 column vector of cell face locations
% Grid.dof = Nx by 1 column vector from 1 to N containig the degrees of freedom, i.e. cell numbers
% Grid.dof_xmin  = scalar cell degree of freedom corrsponding to the left boundary
% Grid.dof_xmax  = scalar cell degree of freedom corrsponding to the right boundary
% Grid.dof_f_xmin = scalar face degree of freedom corrsponding to the left boundary
% Grid.dof_f_xmax = scalar face degree of freedom corrsponding to the right boundary
% + anything else you might find useful
%
% Example call: 
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10; 
% >> Grid = build_grid(Grid);

%% Set up catesian geometry
if ~isfield(Grid,'xmin'); Grid.xmin = 0;  fprintf('Grid.xmin is not defined and has been set to zero.\n');end
if ~isfield(Grid,'xmax'); Grid.xmax = 10; fprintf('Grid.xmax is not defined and has been set to 10.\n'); end
if ~isfield(Grid,'Nx');   Grid.Nx   = 10; fprintf('Grid.Nx is not defined and has been set to 10.\n');end
Grid.Lx = Grid.xmax - Grid.xmin;    % domain length in x
Grid.dx = Grid.Lx/Grid.Nx;        % dx of the gridblocks
% 
% %% Number for fluxes
Grid.Nfx = Grid.Nx + 1;
% 
% % Set up mesh
% % cell centers 'xc' and cell faces 'xf'   
Grid.xc = [Grid.xmin+Grid.dx/2:Grid.dx:Grid.xmax-Grid.dx/2]'; % x-coords of gridblock centers
Grid.xf = [Grid.xmin:Grid.dx:Grid.xmax]'; % x-coords of gridblock faces
% 
% 
% %% Set up dof vectors
Grid.dof = 1:length(Grid.xc);              % cell centered degree of freedom/gridblock number
Grid.dof_f = 1:length(Grid.xf);           % face degree of freedom/face number
% 
% %% Boundary dof's
% % Boundary cells
Grid.dof_xmin = find(mod(Grid.dof,Grid.Nx) == 1);
%Grid.dof_xmin = find(Grid.xc == Grid.xmin + Grid.dx/2);
Grid.dof_xmax = find(mod(Grid.dof,Grid.Nx) == 0);
% % Boundary faces
Grid.dof_f_xmin = find(mod(Grid.dof_f,Grid.Nfx) == 1);
Grid.dof_f_xmax = find(mod(Grid.dof_f,Grid.Nfx) == 0);
end

function [D,G,I]=build_ops(Grid)
% author: Marc Hesse
% date: 09/08/2014
% description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% Output:
% D = Nx by Nx+1 discrete divergence matrix 
% G = Nx+1 by Nx discrete gradient matrix
% I = Nx by Nx identity matrix
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);

% 1) Build sparse Divergence operator
D = 1/Grid.dx.*spdiags([-ones(Grid.Nx,1) ones(Grid.Nx,1)],[0 1],Grid.Nx,Grid.Nfx);
% 2) Obtain sparse Gradient operator in interior from D
G = -D';
% 3) Set natural (homogeneous Neumann) boundary conditions
dof_f_bnd = [find(G(1,:)) find(G(end,:))]; % all dof's on boundary
G([1 end],dof_f_bnd) = 0;

% Identity
I = speye(Grid.Nx);
end