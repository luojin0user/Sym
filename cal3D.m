basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));
% 当运行main后，保存了obj文件以及ICs文件后，可以直接使用这个文件绘图
% 前提是激励不变，如果改变激励，例如改变线圈位置、电流大小等，需要重新运行main函数

ICs_xoy = load("./mat/xoy/ICs.mat", 'ICs').ICs;
obj_xoy = load("./mat/xoy/obj.mat", 'obj').obj;

ICs_zoy = load("./mat/zoy/ICs.mat", 'ICs').ICs;
obj_zoy = load("./mat/zoy/obj.mat", 'obj').obj;

% 取出这两个方向的磁场强度，然后进行分别求解，然后再叠加

points = 400;
x0 = linspace(0,0.28,points);
y0 = 0.12 .* ones(1, points);
z0 = 0.12 .* ones(1, points);

% 注意输入需要是行向量
[Bx, By_xoy] = obj_xoy.cal_Bx_By(ICs_xoy, x0, y0);
[Bz, By_zoy] = obj_zoy.cal_Bx_By(ICs_zoy, z0, y0);

By = By_xoy + By_zoy;

figure; hold on;
plot(x0, Bx, 'b-');  % 蓝
plot(x0, By, 'r--');  % 红
plot(x0, By_xoy, 'g-.'); % 绿
plot(x0, By_zoy, 'g-..'); % 绿

grid on;