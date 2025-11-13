function [x_segments, y_segments, rect_ids] = split_curve_by_rects(x0, y0, rects)
% 输入:
%   x0, y0 : 曲线坐标 (1×N 或 N×1)
%   rects  : 每行一个矩形 [xmin xmax ymin ymax]
% 输出:
%   x_segments, y_segments : cell 数组，每个 cell 对应一段曲线
%   rect_ids : 每段对应的矩形编号

x_segments = {};
y_segments = {};
rect_ids = [];

for k = 1:size(rects,1)
    xmin = rects(k,1); xmax = rects(k,2);
    ymin = rects(k,3); ymax = rects(k,4);
    
    % 判断每个点是否在矩形内部
    in_rect = (x0 >= xmin & x0 <= xmax & y0 >= ymin & y0 <= ymax);
    
    % 找到进入与离开矩形的索引
    diff_in = diff([0, in_rect(:).', 0]);   % 保证是行向量
    start_idx = find(diff_in == 1);
    end_idx   = find(diff_in == -1) - 1;
    
    % 遍历每个连续段
    for i = 1:length(start_idx)
        x_segments{end+1} = x0(start_idx(i):end_idx(i));
        y_segments{end+1} = y0(start_idx(i):end_idx(i));
        rect_ids(end+1) = k;   % 矩形编号
    end
end
end
