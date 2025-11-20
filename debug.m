basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));

% coil = RectangleCoil(0.04, 0.12, 0.10, 0.20, 0.18, 1600, 5);
% coil.set_Rcoil_loc(0.5, 0.5, 0.5);
coil = RectangleCoil(0.04, 100, 200, 98, 198, 1600, 5);
coil.set_Rcoil_loc(0, 0, 0);
fx = coil.factor_x();

tic
fx(0,0,0)
toc