Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 20;
Grid.ymin = 0; Grid.ymax = 1; Grid.Ny = 11;
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
