% --- 0. 清理工作区 ---
clear;
clc;

% --- 1. 定义基本物理常数和几何参数 (根据文档Table 1和Figure 1精确调整) ---
mu0 = 4 * pi * 1e-7; % 真空磁导率 [H/m]

% 几何参数 (cm -> m)
x = [0, 10, 12, 16, 18, 28] * 0.01; % x坐标 [mm]
y = [0, 10, 14, 24] * 0.01; % y坐标 [mm]
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
mu_air = mu0 * mu_r_air;
mu_iron = mu0 * mu_r_iron;

% 区域磁导率 (根据文档定义)
mu1 = mu_air;
mu2 = mu_air;
mu3 = mu_air;
mu4 = mu_air;
mu5 = mu_iron; % Region 5 是铁芯
mu6 = mu_air; % Region 6 是线圈导体区域 (铜线，但考虑为空气磁导率)
mu7 = mu_air; % Region 7 是线圈导体区域 (铜线，但考虑为空气磁导率)

% 谐波阶数 (根据文档Table 1调整)
%H1max = 300; H2max = 300; H3max = 107; H4max = 107; H5max = 43; H6max = 21; H7max = 21;
%N3max = 268; N4max = 268; N5max = 268; N6max = 268; N7max = 268;
H1max = 60; H2max = 60; H3max = 60; H4max = 60; H5max = 60; H6max = 60; H7max = 60;
N3max = 60; N4max = 60; N5max = 60; N6max = 60; N7max = 60;


% 电流密度 (根据文档Table 1调整)
Nt = 1600; I = 5; Sc = 800 * 1e-6; % 匝数，电流，导体表面积 (mm^2 -> m^2)
Jz6 = Nt * I / Sc; % Region 6 的电流密度
Jz7 = -Nt * I / Sc; % Region 7 的电流密度 (反向)

% --- 2. 计算 AzP 和 BxP/ByP (泊松方程的特定解) ---
% 这些值在计算 Q 和 ES 项时，需要在特定 y 坐标处评估。
% AzP6 (27d): - 1/2 * mu6 * Jz6 * y^2
% AzP7 (32d): - 1/2 * mu7 * Jz7 * y^2
% BxP6 (28d): - mu6 * Jz6 * y
% ByP6 (29d): 0
% BxP7 (33d): -mu7 * Jz7 * y
% ByP7 (34d): 0

% 在 y2 和 y3 处评估的特定解值 便于后续代入
global AzP6_y2 AzP6_y3 BxP6_y2 BxP6_y3 ByP6_y2 ByP6_y3
global AzP7_y2 AzP7_y3 BxP7_y2 BxP7_y3 ByP7_y2 ByP7_y3

AzP6_y2 = - 0.5 * mu6 * Jz6 * y2^2;
AzP6_y3 = 0.5 * mu6 * Jz6 * y3^2;
BxP6_y2 = -mu6 * Jz6 * y2;
BxP6_y3 = -mu6 * Jz6 * y3;
ByP6_y2 = 0; % From (29d)
ByP6_y3 = 0; % From (29d)

AzP7_y2 = 0.5 * mu7 * Jz7 * y2^2;
AzP7_y3 = 0.5 * mu7 * Jz7 * y3^2;
BxP7_y2 = -mu7 * Jz7 * y2;
BxP7_y3 = -mu7 * Jz7 * y3;
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
    Bx = 0;
    By = 0;

    switch region_id
        case 1 % Region 1: Laplace's equation, Az imposed on all edges (Eq 4-6)
            % Constants: d1_h1
            for h1_idx = 1:H1max
                h1 = h1_idx;
                d1_h = IC(ic_idx_map.R1_d(h1));
                beta1_h = h1 * pi / Tx1;

                % Az1 (Eq 4)
                Az = Az + d1_h / beta1_h * (sh(beta1_h * (y_coord - y1)) / ch(beta1_h * Ty1)) * sin(beta1_h * (x_coord - x1));

                % Bx1 (Eq 5)
                Bx = Bx + d1_h * (ch(beta1_h * (y_coord - y1)) / ch(beta1_h * Ty1)) * sin(beta1_h * (x_coord - x1));

                % By1 (Eq 6)
                By = By + d1_h * sh(beta1_h * (y_coord - y1)) / ch(beta1_h * Ty1) * cos(beta1_h * (x_coord - x1));
            end
            By = -By;
        case 2 % Region 2: Laplace's equation (Eq 8-10)
            % Constants: c2_h2
            for h2_idx = 1:H2max
                h2 = h2_idx;
                c2_h = IC(ic_idx_map.R2_c(h2));
                beta2_h = h2 * pi / Tx2;

                % Az2 (Eq 8)
                Az = Az + c2_h / beta2_h * (sh(beta2_h * (y4 - y_coord)) / ch(beta2_h * Ty2)) * sin(beta2_h * (x_coord - x1));

                % Bx2 (Eq 9)
                Bx = Bx + c2_h * (ch(beta2_h * (y4 - y_coord)) / ch(beta2_h * Ty2)) * sin(beta2_h * (x_coord - x1));

                % By2 (Eq 10)
                By = By + c2_h * sh(beta2_h * (y4 - y_coord)) / ch(beta2_h * Ty2) * cos(beta2_h * (x_coord - x1));
            end

        case 3 % Region 3: Laplace's equation (Eq 12a-14c)
            % Constants: c3_h3, d3_h3, f3_n3
            % Az3 = Az3_x + Az3_y (12a)
            Az_x = 0; Az_y = 0;
            Bx_x = 0; Bx_y = 0;
            By_x = 0; By_y = 0;

            for h3_idx = 1:H3max
                h3 = h3_idx;
                c3_h = IC(ic_idx_map.R3_c(h3));
                d3_h = IC(ic_idx_map.R3_d(h3));
                beta3_h = h3 * pi / Tx3;

                % Az3_x (12b)
                Az_x = Az_x + (c3_h / beta3_h * sh(beta3_h * (y3 - y_coord)) / sh(beta3_h * Ty3) + ...
                               d3_h / beta3_h * sh(beta3_h * (y_coord - y2)) / sh(beta3_h * Ty3)) * sin(beta3_h * (x_coord - x1));
                % Bx3_x (13b)
                Bx_x = Bx_x + (-c3_h * ch(beta3_h * (y3 - y_coord)) / sh(beta3_h * Ty3) + ...
                               d3_h * ch(beta3_h * (y_coord - y2)) / sh(beta3_h * Ty3)) * sin(beta3_h * (x_coord - x1));
                % By3_x (14b)
                By_x = By_x - (c3_h * (sh(beta3_h * (y3 - y_coord)) / sh(beta3_h * Ty3)) + ...
                               d3_h * (sh(beta3_h * (y_coord - y2)) / sh(beta3_h * Ty3))) * cos(beta3_h * (x_coord - x1));
            end
            for n3_idx = 1:N3max
                n3 = n3_idx;
                f3_n = IC(ic_idx_map.R3_f(n3));
                lambda3_n = n3 * pi / Ty3;

                % Az3_y (12c)
                Az_y = Az_y + f3_n / lambda3_n * (sh(lambda3_n * (x_coord - x1)) / sh(lambda3_n * Tx3)) * sin(lambda3_n * (y_coord - y2));
                % Bx3_y (13c)
                Bx_y = Bx_y + f3_n * (sh(lambda3_n * (x_coord - x1)) / sh(lambda3_n * Tx3)) * cos(lambda3_n * (y_coord - y2));
                % By3_y (14c)
                By_y = By_y - f3_n * (ch(lambda3_n * (x_coord - x1)) / sh(lambda3_n * Tx3)) * sin(lambda3_n * (y_coord - y2));
            end
            Az = Az_x + Az_y;
            Bx = Bx_x + Bx_y;
            By = By_x + By_y;

        case 4 % Region 4: Laplace's equation (Eq 17a-19c)
            % Constants: c4_h4, d4_h4, e4_n4
            Az_x = 0; Az_y = 0;
            Bx_x = 0; Bx_y = 0;
            By_x = 0; By_y = 0;

            for h4_idx = 1:H4max
                h4 = h4_idx;
                c4_h = IC(ic_idx_map.R4_c(h4));
                d4_h = IC(ic_idx_map.R4_d(h4));
                beta4_h = h4 * pi / Tx4;

                % Az4_x (17b)
                Az_x = Az_x + (c4_h * sh(beta4_h * (y3 - y_coord)) / sh(beta4_h * Ty4) + ...
                               d4_h * sh(beta4_h * (y_coord - y2)) / sh(beta4_h * Ty4)) * sin(beta4_h * (x_coord - x5));
                % Bx4_x (18b)
                Bx_x = Bx_x + (c4_h * ch(beta4_h * (y3 - y_coord)) / sh(beta4_h * Ty4) - ...
                               d4_h * ch(beta4_h * (y_coord - y2)) / sh(beta4_h * Ty4)) * sin(beta4_h * (x_coord - x5));
                % By4_x (19b)
                By_x = By_x + (c4_h * beta4_h * sh(beta4_h * (y3 - y_coord)) / sh(beta4_h * Ty4) + ...
                               d4_h * beta4_h * sh(beta4_h * (y_coord - y2)) / sh(beta4_h * Ty4)) * cos(beta4_h * (x_coord - x5));
            end
            for n4_idx = 1:N4max
                n4 = n4_idx;
                e4_n = IC(ic_idx_map.R4_e(n4));
                lambda4_n = n4 * pi / Ty4;

                % Az4_y (17c)
                Az_y = Az_y + e4_n * (sh(lambda4_n * (x6 - x_coord)) / sh(lambda4_n * Tx4)) * sin(lambda4_n * (y_coord - y2));
                % Bx4_y (18c)
                Bx_y = Bx_y - e4_n * (lambda4_n / sh(lambda4_n * Tx4)) * ch(lambda4_n * (x6 - x_coord)) * cos(lambda4_n * (y_coord - y2));
                % By4_y (19c)
                By_y = By_y + e4_n * (ch(lambda4_n * (x6 - x_coord)) / sh(lambda4_n * Tx4)) * sin(lambda4_n * (y_coord - y2));
            end
            Az = Az_x + Az_y;
            Bx = Bx_x + Bx_y;
            By = By_x + By_y;

        case 5 % Region 5: Laplace's equation (Eq 22a-24c)
            % Constants: c5_h5, d5_h5, e5_n5, f5_n5
            Az_x = 0; Az_y = 0;
            Bx_x = 0; Bx_y = 0;
            By_x = 0; By_y = 0;

            for h5_idx = 1:H5max
                h5 = h5_idx;
                c5_h = IC(ic_idx_map.R5_c(h5));
                d5_h = IC(ic_idx_map.R5_d(h5));
                beta5_h = h5 * pi / Tx5;

                % Az5_x (22b)
                Az_x = Az_x + (c5_h * sh(beta5_h * (y3 - y_coord)) / sh(beta5_h * Ty5) + ...
                               d5_h * sh(beta5_h * (y_coord - y2)) / sh(beta5_h * Ty5)) * sin(beta5_h * (x_coord - x3));
                % Bx5_x (23b)
                Bx_x = Bx_x + (c5_h * ch(beta5_h * (y3 - y_coord)) / sh(beta5_h * Ty5) - ...
                               d5_h * ch(beta5_h * (y_coord - y2)) / sh(beta5_h * Ty5)) * sin(beta5_h * (x_coord - x3));
                % By5_x (24b)
                By_x = By_x + (c5_h * beta5_h * sh(beta5_h * (y3 - y_coord)) / sh(beta5_h * Ty5) + ...
                               d5_h * beta5_h * sh(beta5_h * (y_coord - y2)) / sh(beta5_h * Ty5)) * cos(beta5_h * (x_coord - x3));
            end
            for n5_idx = 1:N5max
                n5 = n5_idx;
                e5_n = IC(ic_idx_map.R5_e(n5));
                f5_n = IC(ic_idx_map.R5_f(n5));
                lambda5_n = n5 * pi / Ty5;

                % Az5_y (22c)
                Az_y = Az_y + (e5_n * sh(lambda5_n * (x4 - x_coord)) / sh(lambda5_n * Tx5) + ...
                               f5_n * sh(lambda5_n * (x_coord - x3)) / sh(lambda5_n * Tx5)) * sin(lambda5_n * (y_coord - y2));
                % Bx5_y (23c)
                Bx_y = Bx_y - (e5_n * lambda5_n * ch(lambda5_n * (x4 - x_coord)) / sh(lambda5_n * Tx5) - ...
                               f5_n * lambda5_n * ch(lambda5_n * (x_coord - x3)) / sh(lambda5_n * Tx5)) * cos(lambda5_n * (y_coord - y2));
                % By5_y (24c)
                By_y = By_y + (e5_n * ch(lambda5_n * (x4 - x_coord)) / sh(lambda5_n * Tx5) + ...
                               f5_n * ch(lambda5_n * (x_coord - x3)) / sh(lambda5_n * Tx5)) * sin(lambda5_n * (y_coord - y2));
            end
            Az = Az_x + Az_y;
            Bx = Bx_x + Bx_y;
            By = By_x + By_y;

        case 6 % Region 6: Poisson's equation (Eq 27a-29d)
            % Constants: c6_0, c6_h6, d6_0, d6_h6, e6_n6, f6_n6
            % Az6 = Az6_x + Az6_y + AzP6 (27a)
            Az_x = 0; Az_y = 0;
            Bx_x = 0; Bx_y = 0;
            By_x = 0; By_y = 0;

            % h6 = 0 case
            c6_0 = IC(ic_idx_map.R6_c0);
            d6_0 = IC(ic_idx_map.R6_d0);
            Az_x = Az_x + (y3 - y_coord) * c6_0 + (y_coord - y2) * d6_0;
            Bx_x = Bx_x + (-c6_0 + d6_0); % From (28b) for h=0

            for h6_idx = 1:H6max
                h6 = h6_idx;
                c6_h = IC(ic_idx_map.R6_c(h6));
                d6_h = IC(ic_idx_map.R6_d(h6));
                beta6_h = h6 * pi / Tx6;

                % Az6_x (27b)
                Az_x = Az_x + (c6_h * sh(beta6_h * (y3 - y_coord)) / sh(beta6_h * Ty6) + ...
                               d6_h * sh(beta6_h * (y_coord - y2)) / sh(beta6_h * Ty6)) * cos(beta6_h * (x_coord - x2));
                % Bx6_x (28b)
                Bx_x = Bx_x + (c6_h * ch(beta6_h * (y3 - y_coord)) / sh(beta6_h * Ty6) - ...
                               d6_h * ch(beta6_h * (y_coord - y2)) / sh(beta6_h * Ty6)) * cos(beta6_h * (x_coord - x2));
                % By6_x (29b)
                By_x = By_x - (c6_h * beta6_h * sh(beta6_h * (y3 - y_coord)) / sh(beta6_h * Ty6) + ...
                               d6_h * beta6_h * sh(beta6_h * (y_coord - y2)) / sh(beta6_h * Ty6)) * sin(beta6_h * (x_coord - x2));
            end
            for n6_idx = 1:N6max
                n6 = n6_idx;
                e6_n = IC(ic_idx_map.R6_e(n6));
                f6_n = IC(ic_idx_map.R6_f(n6));
                lambda6_n = n6 * pi / Ty6;

                % Az6_y (27c)
                Az_y = Az_y + (e6_n * ch(lambda6_n * (x_coord - x2)) / sh(lambda6_n * Tx6) + ...
                               f6_n * ch(lambda6_n * (x3 - x_coord)) / sh(lambda6_n * Tx6)) * sin(lambda6_n * (y_coord - y2));
                % Bx6_y (28c)
                Bx_y = Bx_y - (e6_n * lambda6_n * sh(lambda6_n * (x_coord - x2)) / sh(lambda6_n * Tx6) - ...
                               f6_n * lambda6_n * sh(lambda6_n * (x3 - x_coord)) / sh(lambda6_n * Tx6)) * cos(lambda6_n * (y_coord - y2));
                % By6_y (29c)
                By_y = By_y + (e6_n * ch(lambda6_n * (x_coord - x2)) / sh(lambda6_n * Tx6) + ...
                               f6_n * ch(lambda6_n * (x3 - x_coord)) / sh(lambda6_n * Tx6)) * sin(lambda6_n * (y_coord - y2));
            end

            % Add particular solution components
            AzP6_val = 0.5 * mu6 * Jz6 * y_coord^2; % AzP6 (27d)
            BxP6_val = -mu6 * Jz6 * y_coord; % BxP6 (28d)
            ByP6_val = 0; % ByP6 (29d)

            Az = Az_x + Az_y + AzP6_val;
            Bx = Bx_x + Bx_y + BxP6_val;
            By = By_x + By_y + ByP6_val;

        case 7 % Region 7: Poisson's equation (Eq 32a-34d)
            % Constants: c7_0, c7_h7, d7_0, d7_h7, e7_n7, f7_n7
            % Az7 = Az7_x + Az7_y + AzP7 (32a)
            Az_x = 0; Az_y = 0;
            Bx_x = 0; Bx_y = 0;
            By_x = 0; By_y = 0;

            % h7 = 0 case
            c7_0 = IC(ic_idx_map.R7_c0);
            d7_0 = IC(ic_idx_map.R7_d0);
            Az_x = Az_x + (y3 - y_coord) * c7_0 + (y_coord - y2) * d7_0;
            Bx_x = Bx_x + (-c7_0 + d7_0); % From (33b) for h=0

            for h7_idx = 1:H7max
                h7 = h7_idx;
                c7_h = IC(ic_idx_map.R7_c(h7));
                d7_h = IC(ic_idx_map.R7_d(h7));
                beta7_h = h7 * pi / Tx7;

                % Az7_x (32b)
                Az_x = Az_x + (c7_h * sh(beta7_h * (y3 - y_coord)) / sh(beta7_h * Ty7) + ...
                               d7_h * sh(beta7_h * (y_coord - y2)) / sh(beta7_h * Ty7)) * cos(beta7_h * (x_coord - x4));
                % Bx7_x (33b)
                Bx_x = Bx_x + (c7_h * ch(beta7_h * (y3 - y_coord)) / sh(beta7_h * Ty7) - ...
                               d7_h * ch(beta7_h * (y_coord - y2)) / sh(beta7_h * Ty7)) * cos(beta7_h * (x_coord - x4));
                % By7_x (34b)
                By_x = By_x - (c7_h * beta7_h * sh(beta7_h * (y3 - y_coord)) / sh(beta7_h * Ty7) + ...
                               d7_h * beta7_h * sh(beta7_h * (y_coord - y2)) / sh(beta7_h * Ty7)) * sin(beta7_h * (x_coord - x4));
            end
            for n7_idx = 1:N7max
                n7 = n7_idx;
                e7_n = IC(ic_idx_map.R7_e(n7));
                f7_n = IC(ic_idx_map.R7_f(n7));
                lambda7_n = n7 * pi / Ty7;

                % Az7_y (32c)
                Az_y = Az_y + (e7_n * ch(lambda7_n * (x_coord - x4)) / sh(lambda7_n * Tx7) + ...
                               f7_n * ch(lambda7_n * (x5 - x_coord)) / sh(lambda7_n * Tx7)) * sin(lambda7_n * (y_coord - y2));
                % Bx7_y (33c)
                Bx_y = Bx_y - (e7_n * lambda7_n * sh(lambda7_n * (x_coord - x4)) / sh(lambda7_n * Tx7) - ...
                               f7_n * lambda7_n * sh(lambda7_n * (x5 - x_coord)) / sh(lambda7_n * Tx7)) * cos(lambda7_n * (y_coord - y2));
                % By7_y (34c)
                By_y = By_y + (e7_n * ch(lambda7_n * (x_coord - x4)) / sh(lambda7_n * Tx7) + ...
                               f7_n * ch(lambda7_n * (x5 - x_coord)) / sh(lambda7_n * Tx7)) * sin(lambda7_n * (y_coord - y2));
            end

            % Add particular solution components
            AzP7_val = 0.5 * mu7 * Jz7 * y_coord^2; % AzP7 (32d)
            BxP7_val = -mu7 * Jz7 * y_coord; % BxP7 (33d)
            ByP7_val = 0; % ByP7 (34d)

            Az = Az_x + Az_y + AzP7_val;
            Bx = Bx_x + Bx_y + BxP7_val;
            By = By_x + By_y + ByP7_val;

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
num_points = 100;
x_path1 = linspace(x1, x6, num_points);
y_path1 = ((y1 + y2) / 2) * ones(1, num_points);
Bx_path1 = zeros(1, num_points);
By_path1 = zeros(1, num_points);
%disp(y_path1(1:20));

fprintf('\n--- 计算 Path 1 上的磁通密度 ---\n');
for i = 1:num_points
    % Path 1 位于 Region 1
    [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 1, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
end

figure;
plot(x_path1 * 100, Bx_path1, 'b-', 'DisplayName', 'Bx (Subdomain Model)');
hold on;
plot(x_path1 * 100, By_path1, 'r--', 'DisplayName', 'By (Subdomain Model)');
xlabel('Length of Path 1 [cm]');
ylabel('Magnetic Flux Density [T]');
title('Magnetic Flux Density along Path 1');
legend('show');
grid on;
hold off;

% Path 2: y = ((y1 + y4) / 2), x from x1 to x6 (Figure 7)
num_points = 100;
x_path1 = linspace(x1, x6, num_points);
y_path1 = ((y1 + y4) / 2) * ones(1, num_points);
Bx_path1 = zeros(1, num_points);
By_path1 = zeros(1, num_points);

fprintf('\n--- 计算 Path 2 上的磁通密度 ---\n');
for i = 1:num_points
    % Path 2 位于 Region 3,6,5,7,4
    if x_path1(i) > x1 &&  x_path1(i) <= x2
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 3, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    elseif x_path1(i) > x2 &&  x_path1(i) <= x3
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 6, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    elseif x_path1(i) > x3 &&  x_path1(i) <= x4
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 5, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    elseif x_path1(i) > x4 &&  x_path1(i) <= x5
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 7, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    elseif x_path1(i) > x5 &&  x_path1(i) <= x6
        [~, Bx_path1(i), By_path1(i)] = calculate_magnetic_field(x_path1(i), y_path1(i), IC, 4, H1max, H2max, H3max, H4max, H5max, H6max, H7max, N3max, N4max, N5max, N6max, N7max, mu6, mu7, Jz6, Jz7);
    end
end

figure;
plot(x_path1 * 100, Bx_path1, 'b-', 'DisplayName', 'Bx (Subdomain Model)');
hold on;
plot(x_path1 * 100, By_path1, 'r--', 'DisplayName', 'By (Subdomain Model)');
xlabel('Length of Path 1 [cm]');
ylabel('Magnetic Flux Density [T]');
title('Magnetic Flux Density along Path 1');
legend('show');
grid on;
hold off;

% --- 辅助函数定义 (sh, ch, coth, csch) ---
function val = sh(x), val = sinh(x); end
function val = ch(x), val = cosh(x); end
function val = coth(x), val = 1 ./ tanh(x); end
function val = csch(x), val = 1 ./ sinh(x); end

