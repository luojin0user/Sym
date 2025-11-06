C = {[], 1:3, [], 'abc', [], [5 6]}; % 1x6 cell

% 找到非空 cell 的逻辑索引
non_empty_logical = ~cellfun(@isempty, C);

% 得到非空 cell 的索引
non_empty_idx = find(non_empty_logical);

disp(non_empty_logical)
disp(non_empty_idx)