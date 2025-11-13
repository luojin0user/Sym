% 曲线
t = linspace(0,10,100);
x0 = t;
y0 = sin(t);

% 矩形区域
rects = [2 4 -0.5 0.5;
    7 9 -1 1];

[x_seg, y_seg, rect_ids] = split_curve_by_rects(x0, y0, rects);

% 可视化
figure; hold on;
plot(x0, y0, 'k--');  % 原曲线
colors = lines(size(rects,1));

for i = 1:length(x_seg)
    plot(x_seg{i}, y_seg{i}, 'Color', colors(rect_ids(i),:), 'LineWidth',2);
end
axis equal;
