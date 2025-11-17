basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));
% 当运行main后，保存了obj文件以及ICs文件后，可以直接使用这个文件绘图
% 前提是激励不变，如果改变激励，例如改变线圈位置、电流大小等，需要重新运行main函数

ICs = load("./mat/zoy/ICs.mat", 'ICs').ICs;
obj = load("./mat/zoy/obj.mat", 'obj').obj;

points = 400;
x0 = linspace(0,0.28,points);
y0 = 0.12 .* ones(1, points);

% 注意输入需要是行向量
[Bx, By] = obj.cal_Bx_By(ICs, x0, y0);

figure; hold on;
plot(x0, Bx, 'b-', 'Color', [0 0 1]); % 蓝色线表示 Bxc
plot(x0, By, 'r--', 'Color', [1 0 0]); % 红色线表示 Byc
grid on;

%% 3D绘制

num_points_x = 200;  % X 上采样数
num_points_y = 200;  % Y 上采样数

x_vals = linspace(0, 0.28, num_points_x);
y_vals = linspace(0, 0.24, num_points_y);

[Bx_3D, By_3D] = deal(zeros(num_points_y, num_points_x));

parfor iy = 1:num_points_y
    y_tmp = y_vals(iy) * ones(1, num_points_x);
    [Bx_3D(iy,:), By_3D(iy,:)] = obj.cal_Bx_By(ICs, x_vals, y_tmp);
end



% ======== 绘制 3D 面图 ========
figure;
[Xgrid, Ygrid] = meshgrid(x_vals, y_vals);

surf(Xgrid, Ygrid, By_3D, 'EdgeColor', 'none');
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('By [T]');
title('3D Magnetic Flux Density By');
colorbar;
shading interp;
view(45, 35);
grid on;


figure;
[Xgrid, Ygrid] = meshgrid(x_vals, y_vals);

surf(Xgrid, Ygrid, Bx_3D, 'EdgeColor', 'none');
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('By [T]');
title('3D Magnetic Flux Density Bx');
colorbar;
shading interp;
view(45, 35);
grid on;
