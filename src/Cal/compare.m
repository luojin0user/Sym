%% 读取 CSV
data = readmatrix('./comsol18.csv');
x_csv = data(:,1);
y_csv = data(:,2);
z_csv = data(:,3);

%% 读取 MAT
load('./mat/All_B18.mat');  % 得到 Xgrid, Ygrid, B_3D

%% 偏移量
x0 = 0.14;   % 根据你的实际情况调整
y0 = 0.14;

%% CSV 坐标平移
x_csv_shift = x_csv + x0;
y_csv_shift = y_csv + y0;


%% 将 CSV 数据插值到 MAT 网格上
Z_csv_interp = griddata(x_csv_shift, y_csv_shift, z_csv, Xgrid, Ygrid);

%% 求差值
Diff = abs(Z_csv_interp - B_3D) ./ Z_csv_interp;

%% 绘制差值图
figure;
surf(Xgrid, Ygrid, Diff);
shading interp;
colorbar;
title('Difference (CSV interpolated - MAT)');
xlabel('x'); ylabel('y'); zlabel('\Delta B');

