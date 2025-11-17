basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));

tic;
% xoy方向 Bx By
plain1 = AllRegions("xoy");  % 创建对象
plain1.input_calculate_area(0, 0.28, 0, 0.24, 1);
plain1.input_current_region(0.1, 0.12, 0.1, 0.14, 1, 5, 1600);
plain1.input_current_region(0.16, 0.18, 0.1, 0.14, 1, -5, 1600);
plain1.pre_process();
plain1.get_all_regions();  % 调用方法生成并显示所有区域

% zoy方向 Bz By
plain2 = AllRegions("zoy");  % 创建对象
plain2.input_calculate_area(0, 0.28, 0, 0.24, 1);
plain2.input_current_region(0.05, 0.07, 0.1, 0.14, 1, 5, 1600);
plain2.input_current_region(0.21, 0.23, 0.1, 0.14, 1, -5, 1600);
plain2.pre_process();
plain2.get_all_regions();  % 调用方法生成并显示所有区域
toc;


