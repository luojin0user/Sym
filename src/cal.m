% --- 0. 清理工作区 ---
clear;
clc;

% --- 1. 定义基本物理常数和几何参数 (根据文档Table 1和Figure 1精确调整) ---
mu0 = 4 .* pi .* 1e-7; % 真空磁导率 [H/m]

% 几何参数 (cm -> m)
x = [0, 10, 12, 16, 18, 28] .* 0.01; % x坐标 [mm]
y = [0, 10, 14, 24] .* 0.01; % y坐标 [mm]
global x1 x2 x3 x4 x5 x6
global y1 y2 y3 y4
x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); x5 = x(5); x6 = x(6);
y1 = y(1); y2 = y(2); y3 = y(3); y4 = y(4);


% 区域宽度和高度
global Tx1 Ty1 Tx2 Ty2 Tx3 Ty3 Tx4 Ty4 Tx5 Ty5 Tx6 Ty6 Tx7 Ty7
Tx1 = x6 - x1; % Region 1, 2 的 x 宽度
Ty1 = y2 - y1; % Region 1 的 y 高度
Tx2 = Tx1;     % Region 2 的 x 宽度
Ty2 = y4 - y3; % Region 2 的 y 高度
Tx3 = x2 - x1; % Region 3 的 x 宽度
Ty3 = y3 - y2; % Region 3 的 y 高度
Tx4 = x6 - x5; % Region 4 的 x 宽度
Ty4 = y3 - y2; % Region 4 的 y 高度
Tx5 = x4 - x3; % Region 5 的 x 宽度
Ty5 = y3 - y2; % Region 5 的 y 高度
Tx6 = x3 - x2; % Region 6 的 x 宽度
Ty6 = y3 - y2; % Region 6 的 y 高度
Tx7 = x5 - x4; % Region 7 的 x 宽度
Ty7 = y3 - y2; % Region 7 的 y 高度

% 材料磁导率
mu_r_air = 1;
mu_r_iron = 1500;
mu_air = mu0 .* mu_r_air;
mu_iron = mu0 .* mu_r_iron;

% 区域磁导率 (根据文档定义)
mu1 = mu_air;
mu2 = mu_air;
mu3 = mu_air;
mu4 = mu_air;
mu5 = mu_air; % Region 5 是铁芯
mu6 = mu_air; % Region 6 是线圈导体区域 (铜线，但考虑为空气磁导率)
mu7 = mu_air; % Region 7 是线圈导体区域 (铜线，但考虑为空气磁导率)

% 谐波阶数 (根据文档Table 1调整)
%H1max = 300; H2max = 300; H3max = 107; H4max = 107; H5max = 43; H6max = 21; H7max = 21;
%N3max = 268; N4max = 268; N5max = 268; N6max = 268; N7max = 268;
H1max = 60; H2max = 60; H3max = 60; H4max = 60; H5max = 60; H6max = 60; H7max = 60;
N3max = 60; N4max = 60; N5max = 60; N6max = 60; N7max = 60;

% H1max = 63; H2max = 63; H3max = 63; H4max = 63; H5max = 63; H6max = 63; H7max = 63;
% N3max = 62; N4max = 62; N5max = 62; N6max = 62; N7max = 62;
% 电流密度 (根据文档Table 1调整)
Nt = 1600; I = 5; Sc = 800 .* 1e-6; % 匝数，电流，导体表面积 (mm^2 -> m^2)
Jz6 = Nt .* I / Sc; % Region 6 的电流密度
Jz7 = -Nt .* I / Sc; % Region 7 的电流密度 (反向)

% --- 2. 计算 AzP 和 BxP/ByP (泊松方程的特定解) ---
% 这些值在计算 Q 和 ES 项时，需要在特定 y 坐标处评估。
% AzP6 (27d): - 1/2 .* mu6 .* Jz6 .* y^2
% AzP7 (32d): - 1/2 .* mu7 .* Jz7 .* y^2
% BxP6 (28d): - mu6 .* Jz6 .* y
% ByP6 (29d): 0
% BxP7 (33d): -mu7 .* Jz7 .* y
% ByP7 (34d): 0

% 在 y2 和 y3 处评估的特定解值 便于后续代入
global AzP6_y2 AzP6_y3 BxP6_y2 BxP6_y3 ByP6_y2 ByP6_y3
global AzP7_y2 AzP7_y3 BxP7_y2 BxP7_y3 ByP7_y2 ByP7_y3

AzP6_y2 = - 0.5 .* mu6 .* Jz6 .* y2^2;
AzP6_y3 = 0.5 .* mu6 .* Jz6 .* y3^2;
BxP6_y2 = -mu6 .* Jz6 .* y2;
BxP6_y3 = -mu6 .* Jz6 .* y3;
ByP6_y2 = 0; % From (29d)
ByP6_y3 = 0; % From (29d)

AzP7_y2 = 0.5 .* mu7 .* Jz7 .* y2^2;
AzP7_y3 = 0.5 .* mu7 .* Jz7 .* y3^2;
BxP7_y2 = -mu7 .* Jz7 .* y2;
BxP7_y3 = -mu7 .* Jz7 .* y3;
ByP7_y2 = 0; % From (34d)
ByP7_y3 = 0; % From (34d)


% --- 3. 定义积分常数 (IC) 向量的维度和索引映射 ---
% IC = [d1_h1; c2_h2; c3_h3; d3_h3; f3_n3; ... ]
% 按照公式 38a-38h 的顺序构建 IC 向量

len_IC1 = H1max; % d1_h1
len_IC2 = H2max; % c2_h2
len_IC3 = H3max + H3max + N3max; % c3_h3, d3_h3, f3_n3
len_IC4 = H4max + H4max + N4max; % c4_h4, d4_h4, e4_n4
len_IC5 = H5max + H5max + N5max + N5max; % c5_h5, d5_h5, e5_n5, f5_n5
len_IC6 = 1 + H6max + 1 + H6max + N6max + N6max; % c6_0, c6_h6, d6_0, d6_h6, e6_n6, f6_n6
len_IC7 = 1 + H7max + 1 + H7max + N7max + N7max; % c7_0, c7_h7, d7_0, d7_h7, e7_n7, f7_n7

Xmax = len_IC1 + len_IC2 + len_IC3 + len_IC4 + len_IC5 + len_IC6 + len_IC7;
fprintf('len_IC1:%d\n',len_IC1);
fprintf('len_IC2:%d\n',len_IC2);
fprintf('len_IC3:%d\n',len_IC3);
fprintf('len_IC4:%d\n',len_IC4);
fprintf('len_IC5:%d\n',len_IC5);
fprintf('len_IC6:%d\n',len_IC6);
fprintf('len_IC7:%d\n',len_IC7);
% 定义全局索引，方便在函数中引用
global ic_idx_map
ic_idx_map = struct();
current_idx = 0;

% Region 1: d1_h1
ic_idx_map.R1_d = (current_idx + 1):(current_idx + len_IC1);
current_idx = current_idx + len_IC1;

% Region 2: c2_h2
ic_idx_map.R2_c = (current_idx + 1):(current_idx + len_IC2);
current_idx = current_idx + len_IC2;

% Region 3: c3_h3, d3_h3, f3_n3
ic_idx_map.R3_c = (current_idx + 1):(current_idx + H3max);
current_idx = current_idx + H3max;
ic_idx_map.R3_d = (current_idx + 1):(current_idx + H3max);
current_idx = current_idx + H3max;
ic_idx_map.R3_f = (current_idx + 1):(current_idx + N3max);
current_idx = current_idx + N3max;

% Region 4: c4_h4, d4_h4, e4_n4
ic_idx_map.R4_c = (current_idx + 1):(current_idx + H4max);
current_idx = current_idx + H4max;
ic_idx_map.R4_d = (current_idx + 1):(current_idx + H4max);
current_idx = current_idx + H4max;
ic_idx_map.R4_e = (current_idx + 1):(current_idx + N4max);
current_idx = current_idx + N4max;

% Region 5: c5_h5, d5_h5, e5_n5, f5_n5
ic_idx_map.R5_c = (current_idx + 1):(current_idx + H5max);
current_idx = current_idx + H5max;
ic_idx_map.R5_d = (current_idx + 1):(current_idx + H5max);
current_idx = current_idx + H5max;
ic_idx_map.R5_e = (current_idx + 1):(current_idx + N5max);
current_idx = current_idx + N5max;
ic_idx_map.R5_f = (current_idx + 1):(current_idx + N5max);
current_idx = current_idx + N5max;

% Region 6: c6_0, c6_h6, d6_0, d6_h6, e6_n6, f6_n6
ic_idx_map.R6_c0 = current_idx + 1;
current_idx = current_idx + 1;
ic_idx_map.R6_c = (current_idx + 1):(current_idx + H6max);
current_idx = current_idx + H6max;
ic_idx_map.R6_d0 = current_idx + 1;
current_idx = current_idx + 1;
ic_idx_map.R6_d = (current_idx + 1):(current_idx + H6max);
current_idx = current_idx + H6max;
ic_idx_map.R6_e = (current_idx + 1):(current_idx + N6max);
current_idx = current_idx + N6max;
ic_idx_map.R6_f = (current_idx + 1):(current_idx + N6max);
current_idx = current_idx + N6max;

% Region 7: c7_0, c7_h7, d7_0, d7_h7, e7_n7, f7_n7
ic_idx_map.R7_c0 = current_idx + 1;
current_idx = current_idx + 1;
ic_idx_map.R7_c = (current_idx + 1):(current_idx + H7max);
current_idx = current_idx + H7max;
ic_idx_map.R7_d0 = current_idx + 1;
current_idx = current_idx + 1;
ic_idx_map.R7_d = (current_idx + 1):(current_idx + H7max);
current_idx = current_idx + H7max;
ic_idx_map.R7_e = (current_idx + 1):(current_idx + N7max);
current_idx = current_idx + N7max;
ic_idx_map.R7_f = (current_idx + 1):(current_idx + N7max);
current_idx = current_idx + N7max;

tmp = load("IC.mat");
IC = tmp.IC;
disp(IC(1:min(20, Xmax)));
% --- 1. 定义磁场计算函数 ---
function [Az, Bx, By] = calculate_magnetic_field(x_coord, y_coord, IC, region_id, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7)
% 根据给定的坐标 (x_coord, y_coord) 和区域 ID (region_id)
% 计算该点的磁矢量势 Az 和磁通密度 Bx, By。
% IC: 包含所有积分常数的向量。
global ic_idx_map
global Tx1 Ty1 Tx2 Ty2 Tx3 Ty3 Tx4 Ty4 Tx5 Ty5 Tx6 Ty6 Tx7 Ty7
global x1 x2 x3 x4 x5 x6
global y1 y2 y3 y4
Az = 0;

switch region_id
    case 1 % Region 1: Laplace's equation, Az imposed on all edges (Eq 4-6)
        % Constants: d1_h1
        h1 = 1:H1max;
        d1_h = IC(ic_idx_map.R1_d)';
        beta1_h = h1 .* pi / Tx1;
        % Bx1 (Eq 5)
        Bxv = arrayfun(@(d1_h, beta1_h) ...
            d1_h .* (ch(beta1_h .* (y_coord - y1)) ./ ch(beta1_h .* Ty1)) .* ...
            sin(beta1_h .* (x_coord - x1)), ...
            d1_h, beta1_h);
        
        % By1 (Eq 6)
        Byv = arrayfun(@(d1_h, beta1_h) ...
            -d1_h .* (sh(beta1_h .* (y_coord - y1)) ./ ch(beta1_h .* Ty1)) .* ...
            cos(beta1_h .* (x_coord - x1)), ...
            d1_h, beta1_h);
        
        Bx = sum(Bxv);
        By = sum(Byv);
        
    case 2 % Region 2: Laplace's equation (Eq 8-10)
        % Constants: c2_h2
        h2 = 1:H2max;
        c2_h = IC(ic_idx_map.R2_c)';
        beta2_h = h2 .* pi / Tx2;
        
        % Bx2 (Eq 9)
        Bxv = arrayfun(@(c2_h, beta2_h) ...
            c2_h .* (ch(beta2_h .* (y4 - y_coord)) ./ ch(beta2_h .* Ty2)) .* ...
            sin(beta2_h .* (x_coord - x1)), ...
            c2_h, beta2_h);
        
        % By2 (Eq 10)
        Byv = arrayfun(@(c2_h, beta2_h) ...
            c2_h .* (sh(beta2_h .* (y4 - y_coord)) ./ ch(beta2_h .* Ty2)) .* ...
            cos(beta2_h .* (x_coord - x1)), ...
            c2_h, beta2_h);
        
        Bx = sum(Bxv);
        By = sum(Byv);
        
    case 3 % Region 3: Laplace's equation (Eq 12a-14c)
        % Constants: c3_h3, d3_h3, f3_n3
        % Az3 = Az3_x + Az3_y (12a)
        
        h3 = 1:H3max;
        c3_h = IC(ic_idx_map.R3_c)';
        d3_h = IC(ic_idx_map.R3_d)';
        beta3_h = h3 .* pi / Tx3;
        
        % Bx3_x (13b)
        Bx_xv = arrayfun(@(c3_h, d3_h, beta3_h) ...
            (-c3_h .* ch(beta3_h .* (y3 - y_coord)) ./ sh(beta3_h .* Ty3) + ...
            d3_h .* ch(beta3_h .* (y_coord - y2)) ./ sh(beta3_h .* Ty3)) .* ...
            sin(beta3_h .* (x_coord - x1)), ...
            c3_h, d3_h, beta3_h);
        
        % By3_x (14b)
        By_xv = arrayfun(@(c3_h, d3_h, beta3_h) ...
            (c3_h .* (sh(beta3_h .* (y3 - y_coord)) ./ sh(beta3_h .* Ty3)) + ...
            d3_h .* (sh(beta3_h .* (y_coord - y2)) ./ sh(beta3_h .* Ty3))) .* ...
            cos(beta3_h .* (x_coord - x1)), ...
            c3_h, d3_h, beta3_h);
        
        n3 = 1:N3max;
        f3_n = IC(ic_idx_map.R3_f)';
        lambda3_n = n3 .* pi / Ty3;
        % Bx3_y (13c)
        Bx_yv = arrayfun(@(f3_n, lambda3_n) ...
            f3_n .* (sh(lambda3_n .* (x_coord - x1)) ./ sh(lambda3_n .* Tx3)) .* ...
            cos(lambda3_n .* (y_coord - y2)), ...
            f3_n, lambda3_n);
        
        % By3_y (14c)
        By_yv = arrayfun(@(f3_n, lambda3_n) ...
            f3_n .* (ch(lambda3_n .* (x_coord - x1)) ./ sh(lambda3_n .* Tx3)) .* ...
            sin(lambda3_n .* (y_coord - y2)), ...
            f3_n, lambda3_n);
        
        Bx = sum(Bx_xv) + sum(Bx_yv);
        By = -sum(By_xv) - sum(By_yv);
        
    case 4 % Region 4: Laplace's equation (Eq 17a-19c)
        % Constants: c4_h4, d4_h4, e4_n4
        h4 = 1:H4max;
        c4_h = IC(ic_idx_map.R4_c)';
        d4_h = IC(ic_idx_map.R4_d)';
        beta4_h = h4 .* pi / Tx4;
        
        % Bx4_x (18b)
        Bx_xv = arrayfun(@(c4_h, d4_h, beta4_h) ...
            (-c4_h .* ch(beta4_h .* (y3 - y_coord)) ./ sh(beta4_h .* Ty4) + ...
            d4_h .* ch(beta4_h .* (y_coord - y2)) ./ sh(beta4_h .* Ty4)) .* ...
            sin(beta4_h .* (x_coord - x5)), ...
            c4_h, d4_h, beta4_h);
        
        % By4_x (19b)
        By_xv = arrayfun(@(c4_h, d4_h, beta4_h) ...
            (c4_h .* sh(beta4_h .* (y3 - y_coord)) ./ sh(beta4_h .* Ty4) + ...
            d4_h .* sh(beta4_h .* (y_coord - y2)) ./ sh(beta4_h .* Ty4)) .* ...
            cos(beta4_h .* (x_coord - x5)), ...
            c4_h, d4_h, beta4_h);
        
        n4 = 1:N4max;
        e4_n = IC(ic_idx_map.R4_e)';
        lambda4_n = n4 .* pi / Ty4;
        
        % Bx4_y (18c)
        Bx_yv = arrayfun(@(e4_n, lambda4_n) ...
            e4_n .* (sh(lambda4_n .* (x6 - x_coord)) ./ sh(lambda4_n .* Tx4)) .* ...
            cos(lambda4_n .* (y_coord - y2)), ...
            e4_n, lambda4_n);
        
        % By4_y (19c)
        By_yv = arrayfun(@(e4_n, lambda4_n) ...
            e4_n .* (ch(lambda4_n .* (x6 - x_coord)) ./ sh(lambda4_n .* Tx4)) .* ...
            sin(lambda4_n .* (y_coord - y2)), ...
            e4_n, lambda4_n);
        
        Bx = sum(Bx_xv) + sum(Bx_yv);
        By = -sum(By_xv) + sum(By_yv);
        
    case 5 % Region 5: Laplace's equation (Eq 22a-24c)
        % Constants: c5_h5, d5_h5, e5_n5, f5_n5
        
        h5 = 1:H5max;
        c5_h = IC(ic_idx_map.R5_c)';
        d5_h = IC(ic_idx_map.R5_d)';
        beta5_h = h5 .* pi / Tx5;
        
        % Bx5_x (23b)
        Bx_xv = arrayfun(@(c5_h, d5_h, beta5_h) ...
            (-c5_h .* ch(beta5_h .* (y3 - y_coord)) ./ sh(beta5_h .* Ty5) + ...
            d5_h .* ch(beta5_h .* (y_coord - y2)) ./ sh(beta5_h .* Ty5)) .* ...
            sin(beta5_h .* (x_coord - x3)), ...
            c5_h, d5_h, beta5_h);
        
        % By5_x (24b)
        By_xv = arrayfun(@(c5_h, d5_h, beta5_h) ...
            ( c5_h .* sh(beta5_h .* (y3 - y_coord)) ./ sh(beta5_h .* Ty5) + ...
            d5_h .* sh(beta5_h .* (y_coord - y2)) ./ sh(beta5_h .* Ty5)) .* ...
            cos(beta5_h .* (x_coord - x3)), ...
            c5_h, d5_h, beta5_h);
        
        n5 = 1:N5max;
        e5_n = IC(ic_idx_map.R5_e)';
        f5_n = IC(ic_idx_map.R5_f)';
        lambda5_n = n5 .* pi / Ty5;
        
        % Bx5_y (23c)
        Bx_yv = arrayfun(@(e5_n, f5_n, lambda5_n) ...
            (e5_n .* sh(lambda5_n .* (x4 - x_coord)) ./ sh(lambda5_n .* Tx5) + ...
            f5_n .* sh(lambda5_n .* (x_coord - x3)) ./ sh(lambda5_n .* Tx5)) .* ...
            cos(lambda5_n .* (y_coord - y2)), ...
            e5_n, f5_n, lambda5_n);
        
        % By5_y (24c)
        By_yv = arrayfun(@(e5_n, f5_n, lambda5_n) ...
            (-e5_n .* ch(lambda5_n .* (x4 - x_coord)) ./ sh(lambda5_n .* Tx5) + ...
            f5_n .* ch(lambda5_n .* (x_coord - x3)) ./ sh(lambda5_n .* Tx5)) .* ...
            sin(lambda5_n .* (y_coord - y2)), ...
            e5_n, f5_n, lambda5_n);
        
        Bx = sum(Bx_xv) + sum(Bx_yv);
        By = -sum(By_xv) - sum(By_yv);
        
    case 6 % Region 6: Poisson's equation (Eq 27a-29d)
        % Constants: c6_0, c6_h6, d6_0, d6_h6, e6_n6, f6_n6
        % Az6 = Az6_x + Az6_y + AzP6 (27a)
        
        % h6 = 0 case
        c6_0 = IC(ic_idx_map.R6_c0);
        d6_0 = IC(ic_idx_map.R6_d0);
        Bx_x0 = (-c6_0 + d6_0); % From (28b) for h=0
        
        
        h6 = 1:H6max;
        c6_h = IC(ic_idx_map.R6_c)';
        d6_h = IC(ic_idx_map.R6_d)';
        beta6_h = h6 .* pi / Tx6;
        
        % Bx6_x (28b)
        Bx_xv = arrayfun(@(c6_h, d6_h, beta6_h) ...
            (-c6_h .* ch(beta6_h .* (y3 - y_coord)) ./ sh(beta6_h .* Ty6) + ...
            d6_h .* ch(beta6_h .* (y_coord - y2)) ./ sh(beta6_h .* Ty6)) .* ...
            cos(beta6_h .* (x_coord - x2)), ...
            c6_h, d6_h, beta6_h);
        
        % By6_x (29b)
        By_xv = arrayfun(@(c6_h, d6_h, beta6_h) ...
            (c6_h .* sh(beta6_h .* (y3 - y_coord)) ./ sh(beta6_h .* Ty6) + ...
            d6_h .* sh(beta6_h .* (y_coord - y2)) ./ sh(beta6_h .* Ty6)) .* ...
            sin(beta6_h .* (x_coord - x2)), ...
            c6_h, d6_h, beta6_h);
        
        n6 = 1:N6max;
        e6_n = IC(ic_idx_map.R6_e)';
        f6_n = IC(ic_idx_map.R6_f)';
        lambda6_n = n6 .* pi / Ty6;
        
        % Bx6_y (28c)
        Bx_yv = arrayfun(@(e6_n, f6_n, lambda6_n) ...
            (e6_n .* ch(lambda6_n .* (x_coord - x2)) ./ sh(lambda6_n .* Tx6) - ...
            f6_n .* ch(lambda6_n .* (x3 - x_coord)) ./ sh(lambda6_n .* Tx6)) .* ...
            cos(lambda6_n .* (y_coord - y2)), ...
            e6_n, f6_n, lambda6_n);
        
        % By6_y (29c)
        By_yv = arrayfun(@(e6_n, f6_n, lambda6_n) ...
            (e6_n .* sh(lambda6_n .* (x_coord - x2)) ./ sh(lambda6_n .* Tx6) + ...
            f6_n .* sh(lambda6_n .* (x3 - x_coord)) ./ sh(lambda6_n .* Tx6)) .* ...
            sin(lambda6_n .* (y_coord - y2)), ...
            e6_n, f6_n, lambda6_n);
        
        
        % Add particular solution components
        BxP6_val = -mu6 .* Jz6 .* y_coord; % BxP6 (28d)
        ByP6_val = 0; % ByP6 (29d)
        
        
        Bx = sum(Bx_xv) - sum(Bx_yv) + BxP6_val + Bx_x0;
        By = sum(By_xv) + sum(By_yv) + ByP6_val;
        
    case 7 % Region 7: Poisson's equation (Eq 32a-34d)
        % Constants: c7_0, c7_h7, d7_0, d7_h7, e7_n7, f7_n7
        % Az7 = Az7_x + Az7_y + AzP7 (32a)
        
        % h7 = 0 case
        c7_0 = IC(ic_idx_map.R7_c0);
        d7_0 = IC(ic_idx_map.R7_d0);
        Bx_x0 = (-c7_0 + d7_0); % From (33b) for h=0
        
        h7 = 1:H7max;
        c7_h = IC(ic_idx_map.R7_c)';
        d7_h = IC(ic_idx_map.R7_d)';
        beta7_h = h7 .* pi / Tx7;
        
        % Bx7_x (33b)
        Bx_xv = arrayfun(@(c7_h, d7_h, beta7_h) ...
            (-c7_h .* ch(beta7_h .* (y3 - y_coord)) ./ sh(beta7_h .* Ty7) + ...
            d7_h .* ch(beta7_h .* (y_coord - y2)) ./ sh(beta7_h .* Ty7)) .* ...
            cos(beta7_h .* (x_coord - x4)), ...
            c7_h, d7_h, beta7_h);
        
        % By7_x (34b)
        By_xv = arrayfun(@(c7_h, d7_h, beta7_h) ...
            (c7_h .* sh(beta7_h .* (y3 - y_coord)) ./ sh(beta7_h .* Ty7) + ...
            d7_h .* sh(beta7_h .* (y_coord - y2)) ./ sh(beta7_h .* Ty7)) .* ...
            sin(beta7_h .* (x_coord - x4)), ...
            c7_h, d7_h, beta7_h);
        
        n7 = 1:N7max;
        e7_n = IC(ic_idx_map.R7_e)';
        f7_n = IC(ic_idx_map.R7_f)';
        lambda7_n = n7 .* pi / Ty7;
        
        
        % Bx7_y (33c)
        Bx_yv = arrayfun(@(e7_n, f7_n, lambda7_n) ...
            (e7_n .* ch(lambda7_n .* (x_coord - x4)) ./ sh(lambda7_n .* Tx7) - ...
            f7_n .* ch(lambda7_n .* (x5 - x_coord)) ./ sh(lambda7_n .* Tx7)) .* ...
            cos(lambda7_n .* (y_coord - y2)), ...
            e7_n, f7_n, lambda7_n);
        
        % By7_y (34c)
        By_yv = arrayfun(@(e7_n, f7_n, lambda7_n) ...
            (e7_n .* sh(lambda7_n .* (x_coord - x4)) ./ sh(lambda7_n .* Tx7) + ...
            f7_n .* sh(lambda7_n .* (x5 - x_coord)) ./ sh(lambda7_n .* Tx7)) .* ...
            sin(lambda7_n .* (y_coord - y2)), ...
            e7_n, f7_n, lambda7_n);
        
        
        % Add particular solution components
        BxP7_val = -mu7 .* Jz7 .* y_coord; % BxP7 (33d)
        ByP7_val = 0; % ByP7 (34d)
        
        Bx = sum(Bx_xv) - sum(Bx_yv) + BxP7_val + Bx_x0;
        By = sum(By_xv) + sum(By_yv) + ByP7_val;
        
    otherwise
        error('Invalid region_id. Please provide a region from 1 to 7.');
end
end

% --- 2. 示例：在特定点计算磁场 ---

% 定义一个测试点 (x, y)
test_x = (x1 + x6) / 2; % 示例：x方向中间
test_y = (y1 + y4) / 2; % 示例：y方向中间

% 定义测试区域 (需要根据 test_x, test_y 确定，这里仅为示例)
% 假设 test_x, test_y 落在 Region 1
test_region_id_1 = 1;
[Az_val_1, Bx_val_1, By_val_1] = calculate_magnetic_field(test_x, test_y, IC, test_region_id_1, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
fprintf('\n--- 在 Region %d (%f, %f) 处计算磁场 ---\n', test_region_id_1, test_x, test_y);
fprintf('Az: %e Wb/m\n', Az_val_1);
fprintf('Bx: %e T\n', Bx_val_1);
fprintf('By: %e T\n', By_val_1);

% 假设 test_x, test_y 落在 Region 6 (线圈区域)
test_x_R6 = (x2 + x3) / 2;
test_y_R6 = (y2 + y3) / 2;
test_region_id_6 = 6;
[Az_val_6, Bx_val_6, By_val_6] = calculate_magnetic_field(test_x_R6, test_y_R6, IC, test_region_id_6, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
fprintf('\n--- 在 Region %d (%f, %f) 处计算磁场 ---\n', test_region_id_6, test_x_R6, test_y_R6);
fprintf('Az: %e Wb/m\n', Az_val_6);
fprintf('Bx: %e T\n', Bx_val_6);
fprintf('By: %e T\n', By_val_6);

% --- 3. 绘制磁场分布 (示例：沿 Path 1 的 Bx) ---
% Path 1: y = ((y1 + y2) / 2), x from x1 to x6 (Figure 7)
num_points = 1000;
x_path1 = linspace(x1, x6, num_points);
y_path1 = ((y1 + y2) / 2) .* ones(1, num_points);
Bx_path1 = zeros(1, num_points);
By_path1 = zeros(1, num_points);
%disp(y_path1(1:20));

fprintf('\n--- 计算 Path 1 上的磁通密度 ---\n');
for i = 1:num_points
    % Path 1 位于 Region 1
    [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 1, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
end

figure;
plot(x_path1 .* num_points, Bx_path1, 'b-', 'DisplayName', 'Bx (Subdomain Model)');
hold on;
plot(x_path1 .* num_points, By_path1, 'r--', 'DisplayName', 'By (Subdomain Model)');
xlabel('Length of Path 1 [cm]');
ylabel('Magnetic Flux Density [T]');
title('Magnetic Flux Density along Path 1');
legend('show');
grid on;
hold off;

% Path 2: y = ((y1 + y4) / 2), x from x1 to x6 (Figure 7)
num_points = 200;
x_path1 = linspace(x1, x6, num_points);
y_path1 = ((y2+y3)/2) .* ones(1, num_points);
Bx_path1 = zeros(1, num_points);
By_path1 = zeros(1, num_points);

fprintf('\n--- 计算 Path 2 上的磁通密度 ---\n');
for i = 1:num_points
    % Path 2 位于 Region 3,6,5,7,4
    if x_path1(i) > x1 &&  x_path1(i) < x2
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 3, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    elseif x_path1(i) > x2 &&  x_path1(i) < x3
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 6, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    elseif x_path1(i) > x3 &&  x_path1(i) < x4
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 5, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    elseif x_path1(i) > x4 &&  x_path1(i) < x5
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 7, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    elseif x_path1(i) > x5 &&  x_path1(i) < x6
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 4, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    end
end

figure;
plot(x_path1, Bx_path1, 'b-', 'DisplayName', 'Bx (Subdomain Model)');
hold on;
plot(x_path1, By_path1, 'r--', 'DisplayName', 'By (Subdomain Model)');
xlabel('Length of Path 2 [cm]');
ylabel('Magnetic Flux Density [T]');
title('Magnetic Flux Density along Path 1');
legend('show');
grid on;
hold off;


% Path 3: x = ((x1 + x2) / 2), y from y1 to y4 (Figure 7)
num_points = 1000;
x_path1 = ((x3 + x4) / 2) .* ones(1, num_points);
y_path1 = linspace(y1, y4, num_points);
Bx_path1 = zeros(1, num_points);
By_path1 = zeros(1, num_points);

fprintf('\n--- 计算 Path 3 上的磁通密度 ---\n');
for i = 1:num_points
    % Path 3 位于 Region 1 3 2
    if y_path1(i) > y1 &&  y_path1(i) <= y2
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 1, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    elseif y_path1(i) > y2 &&  y_path1(i) <= y3
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 5, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    elseif y_path1(i) > y3 &&  y_path1(i) <= y4
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 2, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    end
end

figure;
plot(y_path1, Bx_path1, 'b-', 'DisplayName', 'Bx (Subdomain Model)');
hold on;
plot(y_path1, By_path1, 'r--', 'DisplayName', 'By (Subdomain Model)');
xlabel('Length of Path 3 [cm]');
ylabel('Magnetic Flux Density [T]');
title('Magnetic Flux Density along Path 3');
legend('show');
grid on;
hold off;


% ======== Path 2: Scan all y between y2 and y3 ========
num_points_x = 200;  % X 上采样数
num_points_y = 150;  % Y 上采样数

x_vals = linspace(x1, x6, num_points_x);
y_vals = linspace(y1, y4, num_points_y);

[Bx_3D, By_3D] = deal(zeros(num_points_y, num_points_x));

fprintf('\n=== 扫描 Path 2 区域上的磁通密度 (3D) ===\n');

for iy = 1:num_points_y
    for ix = 1:num_points_x
        x_tmp = x_vals(ix);
        y_tmp = y_vals(iy);
        
        if y_tmp >= y1 && y_tmp <= y2
            region = 1;
        elseif y_tmp >= y3 && y_tmp <= y4
            region = 2;
        else
            if x_tmp >= x1 && x_tmp <= x2
                region = 3;
            elseif x_tmp > x2 && x_tmp <= x3
                region = 6;
            elseif x_tmp > x3 && x_tmp <= x4
                region = 5;
            elseif x_tmp > x4 && x_tmp <= x5
                region = 7;
            elseif x_tmp > x5 && x_tmp <= x6
                region = 4;
            else
                continue
            end
        end
        [~, Bx_3D(iy,ix), By_3D(iy,ix)] = ...
            calculate_magnetic_field(x_tmp, y_tmp, IC, region, ...
            H1max,H2max,H3max,H4max,H5max,H6max,H7max, ...
            N3max,N4max,N5max,N6max,N7max, mu6,mu7, Jz6,Jz7);
    end
end

% ======== 绘制 3D 面图 ========
figure;
[Xgrid, Ygrid] = meshgrid(x_vals, y_vals);

surf(Xgrid, Ygrid, By_3D, 'EdgeColor', 'none');
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('By [T]');
title('3D Magnetic Flux Density By across Path 2 area');
colorbar;
shading interp;
view(45, 35);
grid on;



% --- 辅助函数定义 (sh, ch, coth, csch) ---
function val = sh(x), val = sinh(x); end
function val = ch(x), val = cosh(x); end
function val = coth(x), val = 1 ./ tanh(x); end
function val = csch(x), val = 1 ./ sinh(x); end

