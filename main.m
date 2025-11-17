basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));

tic;
all_regions = AllRegions();  % 创建对象
all_regions.input_calculate_area(0, 0.28, 0, 0.24, 1);
all_regions.input_current_region(0.1, 0.12, 0.1, 0.14, 1, 5, 1600);
all_regions.input_current_region(0.16, 0.18, 0.1, 0.14, 1, -5, 1600);
all_regions.pre_process();
all_regions.get_all_regions();  % 调用方法生成并显示所有区域
toc;


