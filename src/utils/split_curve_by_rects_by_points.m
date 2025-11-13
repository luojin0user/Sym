function [x_segments, y_segments, rect_ids] = split_curve_by_rects_by_points(x0, y0, rects)
% 输入:
%   x0, y0 : 曲线坐标 (1×N 或 N×1)
%   rects  : 每行一个矩形 [xmin xmax ymin ymax]
% 输出:
%   x_segments, y_segments : cell 数组，每个 cell 对应一段曲线
%   rect_ids : 每段对应的矩形编号
% 输出顺序与原曲线点顺序一致

x_segments = {};
y_segments = {};
rect_ids = [];

N = length(x0);
seg_start = 1;

while seg_start <= N
    % 当前点属于哪些矩形
    in_rect = find(all([x0(seg_start) >= rects(:,1), x0(seg_start) <= rects(:,2), ...
        y0(seg_start) >= rects(:,3), y0(seg_start) <= rects(:,4)], 2));
    
    if isempty(in_rect)
        % 如果当前点不在任何矩形中，跳到下一点
        seg_start = seg_start + 1;
        continue
    end
    
    rect_id = in_rect(1);  % 如果属于多个矩形，选择第一个
    seg_end = seg_start;
    
    % 找连续点仍在该矩形内
    while seg_end+1 <= N && x0(seg_end+1) >= rects(rect_id,1) && x0(seg_end+1) <= rects(rect_id,2) ...
            && y0(seg_end+1) >= rects(rect_id,3) && y0(seg_end+1) <= rects(rect_id,4)
        seg_end = seg_end + 1;
    end
    
    % 保存该段
    x_segments{end+1} = x0(seg_start:seg_end);
    y_segments{end+1} = y0(seg_start:seg_end);
    rect_ids(end+1) = rect_id;
    
    % 下一个段的起点
    seg_start = seg_end + 1;
end
end
