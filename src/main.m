basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));

tic;
all_regions = AllRegions();  % 创建对象
all_regions.get_all_regions();  % 调用方法生成并显示所有区域
toc;