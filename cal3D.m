basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));
% 当运行main后，保存了obj文件以及ICs文件后，可以直接使用这个文件绘图
% 前提是激励不变，如果改变激励，例如改变线圈位置、电流大小等，需要重新运行main函数

num_points_x = 100;  % X 上采样数
num_points_y = 100;  % Y 上采样数

x0 = linspace(0,0.28,num_points_x);
y0 = linspace(0,0.28,num_points_y);
z0 = 0.15 .* ones(1, num_points_y);

B_3D = deal(zeros(num_points_y, num_points_x));

parfor iy = 1:num_points_y
    y_tmp = y0(iy) * ones(1, num_points_x);
    [Bx, By, Bz] = calB(x0', y_tmp', z0');
    
    B_3D(:,iy) = sqrt(Bx.^2 + By.^2 + Bz.^2);
    % B_3D(:,iy) = sqrt(2*(Bx_new.^2 + Bz_xoz_new.^2));
    % fprintf("this is %d\n", iy);
end

figure;
[Xgrid, Ygrid] = meshgrid(x0, y0);

save("./mat/All_B18.mat", 'Xgrid', 'Ygrid', 'B_3D');

surf(Xgrid, Ygrid, B_3D, 'EdgeColor', 'none');
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('By [T]');
title('3D Magnetic Flux Density By');
colorbar;
shading interp;
view(45, 35);

grid on;