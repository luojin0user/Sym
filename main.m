basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));

tic;

coil = RectangleCoil(0.08, 0.18, 28, 56, 0.707e-3, 5);
coil.set_Rcoil_loc(0.14, 0.14, 0.03);
coil.set_calculate_area(0, 0.28, 0, 0.28, 0, 0.24, 4*pi*1e-7);
coil.gen_all_regions();



toc;


