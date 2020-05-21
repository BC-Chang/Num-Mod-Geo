function [D,G,I]=build_ops(Grid)
% author: Bernie Chang
% date: Hopefully soon. (03/07/2020)
% description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% Output:
% D = discrete divergence matrix
% G = discrete gradient matrix
% I = identity matrix
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 4;
% >> Grid.ymin = 0; Grid.ymax = 1; Grid.Ny = 3;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);

% this will help
Nx = Grid.Nx; Ny = Grid.Ny; N = Grid.N;
%Nz = Grid.Nz;
if (Nx>1) && (Ny>1)  % 2D case
    % 1D divergence matrices
    Dx = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx, [0 1], Nx, Nx+1); % 1D div-matrix in x-dir
    Dy = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy, [0 1], Ny, Ny+1); % 1D div-matrix in y-dir
    Ix = speye(Nx); % 1D Nx-identity
    Iy = speye(Ny); % 1D Ny-identity

    % 2D Tensor-product divergence matrices
    Dx = kron(Dx,Iy);  % 2D div-matrix in x-dir
    Dy = kron(Ix,Dy);  % 2D div-matrix in y-dir

    % Complete 2D divergence
    D = [Dx, Dy];
    % Boundary faces
    dof_f_bnd = [Grid.dof_f_xmin;Grid.dof_f_xmax;Grid.dof_f_ymin;Grid.dof_f_ymax];
    
elseif (Nx>1) && (Ny==1)
    D = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1);
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax];   % boundary faces
elseif (Nx==1) && (Ny>1)
    D = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1);
    dof_f_bnd = [Grid.dof_f_ymin; Grid.dof_f_ymax];   % boundary faces
end

%% Gradient
% adjoint relation
G = -D';
% natural bc's
G(dof_f_bnd,:) = 0;

%% Identity
I = speye(N);

R0 = @(d) spdiags(1./Grid.xc.^(d-1),0,Grid.Nx,Grid.Nx);
R1 = @(d) spdiags(Grid.xf.^(d-1),0,Grid.Nfx,Grid.Nfx);
Div = @(d) R0(d)*D*R1(d);
%% Modification of the divergence to account for non-cartesian geometry
if strcmp(Grid.geom,'cylindrical_r')
    fprintf('Operators built for 1D cylindrical geometry\n.')
    % Here you need to add the modification of D for cylindrical coords.
    D = Div(2);
elseif strcmp(Grid.geom,'spherical_r')
    fprintf('Operators built for 1D spherical geometry\n.')
    % Here you need to add the modification of D for spherical coordinates
    D = Div(3);
elseif strcmp(Grid.geom,'cartesian')
    % nothing needs to be done
else
    error('Unknown geometry.')
end