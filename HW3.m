clear all
close all
clc
% Define domain and grid
Grid.xmin = 1.3; 
Grid.xmax = 15.4; 
Grid.Nx   = 13;
Grid = build_grid(Grid);

% Define discrete operators
[D,G,I] = build_ops(Grid);

%% Define BC's
Param.dof_dir   = [];
Param.dof_f_dir = [];
Param.g         = [];
Param.dof_neu   = [Grid.dof_xmin;Grid.dof_xmax];
Param.dof_f_neu = [Grid.dof_f_xmin;Grid.dof_f_xmax];
Param.qb        = [-3;1];

[B,N,fn] = build_bnd(Param,Grid,I);