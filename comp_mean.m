function [Kd] = comp_mean(K,p,kvkh,Grid)
% Description:
% Takes coefficient field, K, defined at the cell centers and computes the
% mean specified by the power, p and returns it in a sparse diagonal
% matrix, Kd.

% Input:
% K = Ny by Nx matrix of cell centered values
% p = power of the generalized mean
%       1 (arithmetic mean)
%      -1 (harmonic mean)
% kvkh = ratio of vertical to horizontal conductivity/permeability (anisotropy)
% Grid = structure containing information about the grid.

% Output:
% Kd = Nf by Nf diagonal matrix of power means at the cell faces.

% Example call:
% K = @(x,y) 1+x.^3+sqrt(y);
% Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 5;
% Grid.ymin = 0; Grid.ymax = 1; Grid.Ny = 4;
% Grid = build_grid(Grid);
% [Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
% Kd = comp_mean(K(Xc,Yc),1,1,Grid);

if (p == -1) | (p == 1)
    if (Grid.Nx == Grid.N) || (Grid.Ny == Grid.N) % 1D
        mean = zeros(Grid.N+1,1);
        if Grid.Ny > Grid.Nx % y-direction: assumed column vector
            if isrow(K); error('K must be column vector for 1D grid in y-direction.'); end
            mean(2:Grid.Ny) = sum(.5*[K(1:Grid.Ny-1),K(2:Grid.Ny)].^p,2).^(1/p);
        else                 % x-direction: assumed row vector
            if iscolumn(K); error('K must be row vector for 1D grid in x-direction.'); end
            mean(2:Grid.Nx) = sum(.5*[K(1:Grid.Nx-1);K(2:Grid.Nx)].^p,1).^(1/p);
        end
        Kd = spdiags([mean],[0],Grid.N+1,Grid.N+1);
    elseif (Grid.N > Grid.Nx) | (Grid.N > Grid.Ny) % 2D        
        mean_x = zeros(Grid.Ny,Grid.Nx+1); % initialize as in 1D case
        mean_x(:,2:Grid.Nx) = (0.5*(K(:,2:Grid.Nx).^p+K(:,1:Grid.Nx-1).^p)).^(1/p);
        
        mean_y = zeros(Grid.Ny+1,Grid.Nx); % initialize as in 1D case
        mean_y(2:Grid.Ny,:) = (0.5*(K(2:Grid.Ny,:).^p+K(1:Grid.Ny-1,:).^p)).^(1/p);
        
        Kd = spdiags([mean_x(:);kvkh*mean_y(:)],0,Grid.Nf,Grid.Nf);
    else
        error('3D permeability is not implemented')
    end
else
    error('This power does not have significance.')
end
end



