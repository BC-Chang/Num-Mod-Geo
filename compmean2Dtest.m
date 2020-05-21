K = @(x,y) 1+x.^3+sqrt(y);
Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 5;
Grid.ymin = 0; Grid.ymax = 1; Grid.Ny = 4;
Grid = build_grid(Grid);
[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
Kd = comp_mean(K(Xc,Yc),-1,1,Grid);