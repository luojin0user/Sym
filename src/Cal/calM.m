function M = calM(coilConfig, dx, dz, Ip)
% coilConfig:
%   outer = [xL, yL, zL] 外尺寸
%   inner = [xS, yS, zS] 内尺寸
%   Nturn = 每层匝数
%   NzLayer = z方向层数
%   dx,dz = (x,y),z方向的dx，dx，dy相同
%
% 返回互感 M

%% 1. 生成目标线圈离散点（中心对称中空矩形）
[xGrid, yGrid, zGrid, dx_coil] = generateCoilPointsByNum(coilConfig, dx, dz);
grid_size = size(xGrid);
% 展平成列向量方便计算
xVec = xGrid(:);
yVec = yGrid(:);
zVec = zGrid(:);

%% 2. 计算磁场（批量向量化）
[~, ~, Bz] = calB(xVec, yVec, zVec); % 返回列向量

%% 3. 准备每匝对应的点索引
% BxGrid = reshape(Bx, grid_size);
% ByGrid = reshape(By, grid_size);
BzGrid = reshape(Bz, grid_size);

% 对每一匝线圈取到对应这个位置的磁感应强度
dA = dx * dx;
sum_phi = 0;

dz_num = round(coilConfig.thick / dz);
Nt_per_dz = coilConfig.layer / dz_num;
Nt_per_dx = coilConfig.Nt_per_layer / dx_coil;
for i=1:dz_num
    % z方向
    sum_i = 0;
    for j=1:dx_coil
        % x,y方向
        xe = grid_size(1) - j;
        ye = grid_size(2) - j;
        Bzn = BzGrid(j+1:xe,j+1:ye,i);
        
        dphi = Bzn .* dA;
        sum_i = sum_i + sum(dphi, "all") * Nt_per_dx;
    end
    sum_phi = sum_phi + sum_i * Nt_per_dz;
end

M = sum_phi / Ip;

end

%% 辅助函数：按点数生成离散网格
function [xGrid, yGrid, zGrid, dx_coil] = generateCoilPointsByNum(coil, dx, dz)

% 首先，需要知道线圈分割为多少网格
width = coil.x_cl2 - coil.x_cl1;    % 线圈的宽度
dx_coil = round(width / dx);

% x方向分布
x1 = linspace(coil.x_cl1, coil.x_cl2, dx_coil);
x2 = linspace(coil.x_cr2, coil.x_cr1, dx_coil);

% y方向分布
y1 = linspace(coil.y_cl1, coil.y_cl2, dx_coil);
y2 = linspace(coil.y_cr2, coil.y_cr1, dx_coil);

% 然后需要知道中间空间被分为多少网格
% 这里是从外圈到内圈
d_air_x = x2-x1;
d_air_y = y2-y1;

dx_air_x_num = round(d_air_x / dx);
dx_air_y_num = round(d_air_y / dx);

% 所有的分区
dx_all = cell(dx_coil, 1);
dy_all = cell(dx_coil, 1);

for i=1:dx_coil
    dx_all{i} = linspace(x1(i), x2(i), dx_air_x_num(i));
    dy_all{i} = linspace(y1(i), y2(i), dx_air_y_num(i));
end

% z方向分层
dz_num = round(coil.thick / dz);
zVec = linspace(coil.z_b, coil.z_t, dz_num);


% 生成网格
[X, Y, Z] = meshgrid(dx_all{1}, dy_all{1}, zVec);

xGrid = X; yGrid = Y; zGrid = Z;
end
