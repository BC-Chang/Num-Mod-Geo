function [PSI,psi_min,psi_max] = comp_streamfun(q,Grid)
% author: Marc Hesse
% date: 30 Oct 2014, 12 June 2015, 9 July 2015
% Description: The streamfunction of a numerical solution is computed at 
%              the cell corners given the vector of numerical fluxes across 
%              the cell faces.
% Input: q = [qx;qy] Nf by 1 vector of all face fluxes. x-fluxes first.
%        Grid = structure containing all information about the grid.
% Output: PSI = Ny+1 by Nx+1 matrix containing the values of the
%                     streamfunction along the cell corners.
%         psi_min = minimum of PSI
%         psi_max = maximum of PSI


% Integrate horizontally along ymin boundary first then vertically into the domain.
Qymin = Grid.dx*[0 cumsum(q(Grid.dof_f_ymin))'];  % Integral of flow into ymin boundary, this is a vector
Qx = Grid.dy*reshape(q(1:Grid.Nfx),Grid.Ny,Grid.Nx+1); % Horizontal fluxes in Ny by Nx+1 matrix
PSI = cumsum([-Qymin; Qx]); % integrals into domain with Qymin as initial value
psi_min = min(PSI(:));  psi_max = max(PSI(:));
