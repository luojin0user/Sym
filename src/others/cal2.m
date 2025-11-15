basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));

regions = RegionsInput();
regions.set_calculate_area(0, 0.28, 0, 0.24, 1);
regions.set_current_regions(0.1, 0.12, 0.1, 0.14, 1, 1e7);
regions.set_current_regions(0.16, 0.18, 0.1, 0.14, 1, -1e7);

regions.divide_regions();

regions.findNeighbors();

regions.polt_divided_regions();