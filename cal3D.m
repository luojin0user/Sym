basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));
% 当运行main后，保存了obj文件以及ICs文件后，可以直接使用这个文件绘图
% 前提是激励不变，如果改变激励，例如改变线圈位置、电流大小等，需要重新运行main函数

num_points_x = 100;  % X 上采样数
num_points_y = 100;  % Y 上采样数

x0 = linspace(0,0.28,num_points_x);
y0 = linspace(0,0.28,num_points_y);
z0 = 0.15;

% 原始 x0, y0
[xGrid, yGrid] = meshgrid(x0, y0);

% 展开成列向量 (num_points_x * num_points_y) x 1
xv = xGrid(:);
yv = yGrid(:);
zv = z0 * ones(size(xv));    % 统一展开为列向量

% 一次性计算所有点 —— calB 只接收列向量，所以完全符合要求
[Bxv, Byv, Bzv] = calB(xv, yv, zv);

% reshape 回二维矩阵 (对应 y × x)
Bx = reshape(Bxv, size(xGrid));
By = reshape(Byv, size(xGrid));
Bz = reshape(Bzv, size(xGrid));

% 磁场大小
B_3D = sqrt(Bx.^2 + By.^2 + Bz.^2);


figure;

save("./mat/All_B18.mat", 'xGrid', 'yGrid', 'B_3D');

surf(xGrid, yGrid, B_3D, 'EdgeColor', 'none');
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('By [T]');
title('3D Magnetic Flux Density By');
colorbar;
shading interp;
view(45, 35);

grid on;