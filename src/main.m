addpath(genpath(fullfile(pwd, '..\src')));
addpath(genpath(fullfile(pwd, '..\src\region\')));
addpath(genpath(fullfile(pwd, '..\src\boundary\')));

run('global_val.m')


all_regions = AllRegions();  % 创建对象
all_regions.get_all_regions();  % 调用方法生成并显示所有区域
