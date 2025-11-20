basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));
% 当运行main后，保存了obj文件以及ICs文件后，可以直接使用这个文件绘图
% 前提是激励不变，如果改变激励，例如改变线圈位置、电流大小等，需要重新运行main函数

ICs_xoz = load("./mat/xoz/ICs.mat", 'ICs').ICs;
obj_xoz = load("./mat/xoz/obj.mat", 'obj').obj;

ICs_yoz = load("./mat/yoz/ICs.mat", 'ICs').ICs;
obj_yoz = load("./mat/yoz/obj.mat", 'obj').obj;

coil = load("./mat/RCoil.mat", 'obj').obj;

% 取出这两个方向的磁场强度，然后进行分别求解，然后再叠加
%{
points = 1000;
x0 = linspace(0,0.28,points);
y0 = 0.12 .* ones(1, points);
z0 = 0.12 .* ones(1, points);

% 注意输入需要是行向量
[Bx, Bz_xoz] = obj_xoz.cal_Bx_By(ICs_xoz, x0, y0);
[By, Bz_yoz] = obj_yoz.cal_Bx_By(ICs_yoz, z0, y0);
Bz = Bz_xoz + Bz_yoz;
B_3D = sqrt(Bx.^2 + By.^2 + Bz.^2);

figure; hold on;
plot(x0, B_3D, 'b-');  % 蓝
%}
% figure; hold on;
% plot(x0, Bx, 'b-');  % 蓝
% plot(x0, By, 'r--');  % 红
% plot(x0, By_xoy, 'g-.'); % 绿
% plot(x0, By_zoy, 'g-..'); % 绿

factor_x = coil.factor_x();
factor_y = coil.factor_y();

num_points_x = 100;  % X 上采样数
num_points_y = 100;  % Y 上采样数

x0 = linspace(0,0.28,num_points_x);
y0 = linspace(0,0.4,num_points_y);
z0 = 0.12 .* ones(1, num_points_y);

B_3D = deal(zeros(num_points_y, num_points_x));

parfor iy = 1:num_points_y
    y_tmp = y0(iy) * ones(1, num_points_x);
    [Bx, Bz_xoz] = obj_xoz.cal_Bx_By(ICs_xoz, x0, z0);
    [By, Bz_yoz] = obj_yoz.cal_Bx_By(ICs_yoz, y_tmp, z0);
    
    fx = factor_x(x0', y_tmp', z0');
    fy = factor_y(x0', y_tmp', z0');
    Bx_new = Bx .* fx;
    Bz_xoz_new = Bz_xoz .* fx;

    By_new = By .* fy;
    Bz_yoz_new = Bz_yoz .* fy;

    Bz_new = Bz_xoz_new + Bz_yoz_new;

    % Bz = Bz_xoz + Bz_yoz;
    % B_3D(:,iy) = sqrt(Bx.^2 + By.^2 + Bz.^2);
    B_3D(:,iy) = sqrt(2*(Bx_new.^2 + By_new.^2 + Bz_new.^2));
    % B_3D(:,iy) = Bx_new + By_new + Bz_new;
    % fprintf("this is %d\n", iy);
end

figure;
[Xgrid, Ygrid] = meshgrid(x0, y0);

surf(Xgrid, Ygrid, B_3D, 'EdgeColor', 'none');
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('By [T]');
title('3D Magnetic Flux Density By');
colorbar;
shading interp;
view(45, 35);

grid on;