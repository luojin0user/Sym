basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));
figure;
ax = axes;
hold(ax, "on");

xlim(ax, [0 0.28]);
ylim(ax, [0 0.28]);
zlim(ax, [0 0.24]);
grid(ax, "on");
ax.Box = "on";

coil2 = RectangleCoil(0.08, 0.18, 28, 56, 0.707e-3, 5);
coil2.set_Rcoil_loc(0.14, 0.14, 0.15);
coil2.plot3D(ax);


TXCoil = load("./mat/RCoil.mat").obj;
TXCoil.plot3D(ax);

M = calM(coil2, 2.828e-3, 2.828e-3);
disp(M);
