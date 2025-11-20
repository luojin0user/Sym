classdef RegionsInput < handle
    properties
        special_region_area    % 定义特殊区域的位置
        special_region_property    % 定义特殊区域的性质，例如mu_0 I_r N_t
        air_mu_r                % 其他区域的mu_r
        
        all_area           % 所有需要计算的区域坐标
        
        divided_rects       % 划分好的区域
        all_casetype        % 所有区域的类型
        
        % 这里是所有区域的边界，每一个都是一个1xn的cell，每一列都是一个数组储存对应区域的邻接区域
        all_lefts
        all_rights
        all_tops
        all_bottoms

        % 所有的边界种类
        all_BC_types

        % 区域数量
        regions_num
        special_regions_num

        H_max = 60
        N_max = 60

        all_H_max
        all_N_max
        all_mu_r
        all_J_r
    end
    
    methods
        function obj = RegionsInput()
            % 空构造函数
        end
        
        % 这里输入特殊的区域，例如线圈区域、铁磁区域或者铝区域等，需要输入对应的坐标以及J_r
        % 设置电流区域
        function set_current_regions(obj, xl, xr, yb, yt, mu_r, J_r)
            row = [xl, xr, yb, yt];
            obj.special_region_area(end+1, :) = row;
            
            row = [mu_r, J_r];
            obj.special_region_property(end+1, :) = row;
        end
        
        % 设置需要计算的区域，也就是全体区域
        function set_calculate_area(obj, xl, xr, yb, yt, mu_r)
            % 输入的xl,xr,yb,yt分别是左侧x坐标，右侧x坐标，下侧y坐标，上侧y坐标，mu_r指的是这个区域的相对磁导率，一般为1
            row = [xl, xr, yb, yt];
            obj.all_area = row;
            
            obj.air_mu_r = mu_r;
        end
        
        % 分割区域，用于送给后续进行计算
        function divide_regions(obj)
            % S: 1x4  [xL xR yB yT]
            % subregions: cell array, each 1x4

            S = obj.all_area;
            subregions = obj.special_region_area;
            obj.special_regions_num = size(subregions, 1);
            
            % 验证输入
            if isempty(subregions)
                rects = [S(1) S(2) S(3) S(4)];
                return;
            end
            
            % 1) 生成 y 边界（S 的上下 + 每个子区域的上下）
            y_edges = [S(3); S(4); subregions(:,3); subregions(:,4)];
            y_edges = unique(y_edges);
            y_edges = sort(y_edges);
            
            % 2) 准备所有可能的 x 边界（用于中间需要分割的带）
            all_x_edges = unique([S(1); S(2); subregions(:,1); subregions(:,2)]);
            all_x_edges = sort(all_x_edges);
            
            rects = [];
            
            % 3) 对每个 y 带处理
            for j = 1:length(y_edges)-1
                yB = y_edges(j);
                yT = y_edges(j+1);
                
                % 检查是否有子区域与该带交叠（只要有任意交叠即视作被穿过）
                active = 0;
                for k = 1:size(subregions,1)
                    r = subregions(k,:);
                    if (r(3) == yB) && (r(4) == yT)   % 有垂直交叠
                        active = true;
                        break;
                    end
                end
                
                if ~active
                    % 没有子区域穿过：整条作为一个矩形（上下边界对齐）
                    rects(end+1, :) = [S(1), S(2), yB, yT];
                else
                    % 有子区域穿过：按所有 x 边界分割（避免左右合并）
                    for i = 1:length(all_x_edges)-1
                        xL = all_x_edges(i);
                        xR = all_x_edges(i+1);
                        % 仅保留在 S 内的区段
                        if xL >= S(1) && xR <= S(2)
                            % 当前区域如果是有源区域则需要特殊处理
                            s = false;
                            for k = 1:size(subregions,1)
                                r = subregions(k,:);
                                if (r(1) == xL) && (r(2) == xR)
                                    s = true;
                                    break;
                                end
                            end
                            if ~s
                                rects(end+1, :) = [xL, xR, yB, yT];
                            else
                                
                            end
                        end
                    end
                end
            end
            
            
            % 将有源子区域放在最后
            rects = [rects; subregions];
            
            obj.divided_rects = rects;
            obj.regions_num = size(rects, 1);
        end
        
        % 绘制分割好的区域
        function polt_divided_regions(obj)
            figure; hold on; axis equal;
            rects = obj.divided_rects;
            
            % 绘制切分后的矩形（细边框）
            for k = 1:size(rects,1)
                xL = rects(k,1); xR = rects(k,2);
                yB = rects(k,3); yT = rects(k,4);
                rectangle('Position',[xL yB xR-xL yT-yB],'EdgeColor',[0.2 0.2 0.2]);
            end
            
            subregions = obj.special_region_area;
            
            % 如果提供子区域，则用填充高亮显示（半透明）
            if nargin > 1 && ~isempty(subregions)
                cmap = lines(size(subregions,1));
                for k = 1:size(subregions,1)
                    r = subregions(k,:);
                    % 使用 patch 绘制半透明矩形
                    px = [r(1) r(2) r(2) r(1)];
                    py = [r(3) r(3) r(4) r(4)];
                    p = patch(px, py, cmap(mod(k-1,size(cmap,1))+1,:), 'FaceAlpha', 0.25, 'EdgeColor','none');
                end
            end
            
            xlabel('X'); ylabel('Y'); title('Subdivision');
            hold off;
        end
        
        
        function findNeighbors(obj)
            % regions: N x 4 数组，每行：[xL xR yB yT]
            regions = obj.divided_rects;
            current_regions = obj.special_region_area;

            N = size(regions,1);
            i_current = N - size(current_regions, 1) + 1; % 由此之后就是电流区域
            
            casetypes = cell(N, 1);

            left_neighbors   = cell(N,1);
            right_neighbors  = cell(N,1);
            top_neighbors    = cell(N,1);
            bottom_neighbors = cell(N,1);

            regions_bc_type = cell(N,1);
            
            for i = 1:N

            % 初始设置为真正空数组
            left_neighbors{i}   = [];
            right_neighbors{i}  = [];
            top_neighbors{i}    = [];
            bottom_neighbors{i} = [];
        
                xiL = regions(i,1); xiR = regions(i,2);
                yiB = regions(i,3); yiT = regions(i,4);
                
                % 遍历其它区域
                for j = 1:N
                    if i == j, continue; end
                    
                    xjL = regions(j,1); xjR = regions(j,2);
                    yjB = regions(j,3); yjT = regions(j,4);
                    
                    % ============================
                    % 左邻区域：区域 j 的右边界 = 区域 i 的左边界
                    % 并且 y 区间重叠
                    % ============================
                    if abs(xjR - xiL) < 1e-12 && obj.intervalsOverlap([yiB yiT], [yjB yjT])
                        left_neighbors{i}(end+1) = j;
                    end
                    
                    % ============================
                    % 右邻区域：区域 j 的左边界 = 区域 i 的右边界
                    % ============================
                    if abs(xjL - xiR) < 1e-12 && obj.intervalsOverlap([yiB yiT], [yjB yjT])
                        right_neighbors{i}(end+1) = j;
                    end
                    
                    % ============================
                    % 下邻区域：区域 j 的上边界 = 区域 i 的下边界
                    % ============================
                    if abs(yjT - yiB) < 1e-12 && obj.intervalsOverlap([xiL xiR], [xjL xjR])
                        bottom_neighbors{i}(end+1) = j;
                    end
                    
                    % ============================
                    % 上邻区域：区域 j 的下边界 = 区域 i 的上边界
                    % ============================
                    if abs(yjB - yiT) < 1e-12 && obj.intervalsOverlap([xiL xiR], [xjL xjR])
                        top_neighbors{i}(end+1) = j;
                    end

                    % 处理这个区域种类，目前只考虑这三类
                    if length(top_neighbors{i}) > 1 || length(bottom_neighbors{i}) > 1
                        regions_bc_type{i} = BC_TYPE.BBAA;

                        if isempty(top_neighbors{i}) || isempty(bottom_neighbors{i})
                            % 如果上下有一个大于1，而且还有一个为0，那么区域种类一定是BTAir
                            casetypes{i} = CaseType.BTAir;
                        else
                            % 如果不是，那么就是AllayAir
                            casetypes{i} = CaseType.AlleyAir;
                        end

                    elseif i >= i_current
                        regions_bc_type{i} = BC_TYPE.AABB;
                        casetypes{i} = CaseType.FerriteCurrent;
                    else
                        regions_bc_type{i} = BC_TYPE.AAAA;
                        casetypes{i} = CaseType.NormalAir;
                    end
                end
            end

            obj.all_lefts = left_neighbors;
            obj.all_rights = right_neighbors;
            obj.all_tops = top_neighbors;
            obj.all_bottoms = bottom_neighbors;

            obj.all_BC_types = regions_bc_type;
            obj.all_casetype = casetypes;
        end
        
        
        function flag = intervalsOverlap(obj, a, b)
            % a = [a1 a2]
            % b = [b1 b2]
            flag = ~(a(2) <= b(1) || b(2) <= a(1));
        end
        
        % 计算其他所有信息，例如H_max等
        function cal_other_info(obj)
            r_num = obj.regions_num;
            sr_num = obj.special_regions_num;

            obj.all_H_max = obj.H_max * ones(r_num, 1);
            obj.all_N_max = obj.N_max * ones(r_num, 1);
            
            other_mu_r = obj.air_mu_r * ones(r_num - sr_num, 1);
            other_J_r = zeros(r_num - sr_num, 1);

            special_mu_r = zeros(sr_num, 1);
            special_J_r = zeros(sr_num, 1);

            for i=1:sr_num
                special_mu_r(i) = obj.special_region_property(i,1);
                special_J_r(i) = obj.special_region_property(i,2);
            end

            obj.all_mu_r = [other_mu_r; special_mu_r];
            obj.all_J_r = [other_J_r; special_J_r];

        end

        function [H_max, N_max] = rtn_HN_max(obj)
            H_max = obj.all_H_max;
            N_max = obj.all_N_max;
        end

        function [mu_r, J_r] = rtn_mu_J(obj)
            mu_r = obj.all_mu_r;
            J_r = obj.all_J_r;
        end

        function [regions, region_num, current_regions] = rtn_regions(obj)
            regions = obj.divided_rects;
            region_num = obj.regions_num;
            current_regions = obj.special_region_area;

        end

        function [bctype, casetype] = rtn_types(obj)
            bctype = obj.all_BC_types;
            casetype = obj.all_casetype;
        end

        function [lefts, rights, tops, bottoms] = rtn_boundarys(obj)
            lefts = obj.all_lefts;
            rights = obj.all_rights;
            tops = obj.all_tops;
            bottoms = obj.all_bottoms;
        end

        function current_idx = rtn_current_idx(obj)
            current_idx = (obj.regions_num - obj.special_regions_num + 1):(obj.regions_num);
        end
    end
end