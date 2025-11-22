classdef AllRegions < handle
    properties
        regions  % cell array 存储所有区域
        divide_xoz   % 用于处理输入区域等
        divide_yoz
        
        region_num_xoz % 区域数量
        region_num_yoz % 区域数量
        all_region_num
        all_H_max
        all_N_max
        
        all_mu_r
        all_J_r
        
        len_IC_xoz
        len_IC_yoz
        
        % 区域，每一行是一个矩形区域
        regions_area
        all_BC_types
        all_casetypes
        % 有源区域
        current_regions
        current_regions_idx
        
        % 所有区域的上下左右边界
        all_lefts
        all_rights
        all_tops
        all_bottoms
        
        % 当前平面，xoy或者zoy字符串
        this_plain
    end
    
    methods
        function obj = AllRegions()
            obj.divide_xoz = RegionsInput();
            obj.divide_yoz = RegionsInput();
        end
        
        function get_all_regions(obj)
            % obj.regions = cell(obj.region_num, 1);
            obj.set_all_regions();
            disp("区域内方程计算完毕");
            obj.cal_IC_lens();
            % 在所有的区域内部方程完成运算之后，再调用函数进行所有边界的运算
            
            % 计算所有的边界方程
            [BC_funcs, BC_loc, ES_funcs] = obj.cal_all_BCs();
            disp("计算系数方程完成");
            % 注意这里BCs的顺序
            % 每个BCs{i}代表对应区域i的所有边界方程
            % BCs{i}{j,k}代表当前区域按顺序的{j,k}个方程，例如{2,3}表示区域的第二个临界区域的第三个方程（含有e的方程）
            % 区域内的方程的每一行代表一个邻接区域，按“上下左右”的顺序排列，如上边界有多个就是“上上上上下左右”
            % 每一行的每个元素代表这个邻接区域的方程的c0 c d0 d e f分量，如果为{}则没有该分量，直接跳过
            
            % 拼接ES矩阵
            [ES_xoz, ES_yoz] = obj.splice_ES(ES_funcs);
            disp("计算ES方程完成")
            % save("./mat/ES.mat", 'ES');
            
            % 开始拼接所有的矩阵
            [BC_xoz, BC_yoz] = obj.splice_BC(BC_funcs, BC_loc);
            disp("计算BC方程完成");
            % save("./mat/BC.mat", 'BC');
            
            % IC = BC \ ES;
            IC_xoz = lsqr(BC_xoz, ES_xoz, 1e-6, 1000);
            IC_yoz = lsqr(BC_yoz, ES_yoz, 1e-6, 1000);
            disp("计算IC方程完成");
            % filepath = fullfile("mat", obj.this_plain, "IC.mat");
            % save(filepath, "IC");
            
            [ICs_xoz, ICs_yoz] = obj.split_IC(IC_xoz, IC_yoz);
            save("./mat/ICs", 'ICs_xoz', 'ICs_yoz');
            
            save("./mat/all_regions.mat", 'obj');
            
        end
        
        
        function set_all_regions(obj)
            all_regions = cell(obj.all_region_num, 1);
            parfor i=1:obj.all_region_num
                all_regions{i} = Region(i, obj.all_casetypes{i}, obj.regions_area(i,:), obj.all_BC_types{i}, ...
                    obj.all_tops{i}, obj.all_bottoms{i}, obj.all_lefts{i}, obj.all_rights{i}, obj.current_regions_idx, ...
                    obj.all_H_max(i), obj.all_N_max(i), obj.all_mu_r(i), obj.all_J_r(i), obj.all_region_num);
                all_regions{i}.get_region_solution_func();
            end
            obj.regions = all_regions;
        end
        
        function [BC_funcs, BC_loc, ES_funcs] = cal_all_BCs(obj)
            BC_funcs = cell(obj.all_region_num, 1);
            BC_loc = cell(obj.all_region_num, 1);
            ES_funcs = cell(obj.all_region_num, 1);
            
            parfor i=1:obj.all_region_num    % 并行计算所有边界函数
                [BC_funcs{i}, BC_loc{i}, ES_funcs{i}] = obj.regions{i}.gen_region_coefficient_func(obj.regions);
            end
        end
        
        % 拼接所有的BC
        function [BC_xoz, BC_yoz] = splice_BC(obj, BC_funcs, BC_loc)
            
            N = obj.all_region_num;
            BC_blocks = cell(N, N);
            
            sprintf("开始计算BC矩阵，共%d个\n",N);
            parfor i = 1:N
                for j = 1:N
                    if i == j
                        BC_blocks{i,j} = [];    % 暂时置为空
                        continue;
                    end
                    
                    if obj.regions{i}.all_edge_regions(j)   % 有边界相邻
                        BC_blocks{i,j} = obj.splice_BCxx(i, j, BC_funcs, BC_loc);
                        fprintf("计算完成区域(%d,%d)\n",i,j);
                    else
                        BC_blocks{i,j} = []; % 或保持空
                    end
                end
            end
            
            % 然后分割BC_block
            k = obj.region_num_xoz;
            BC_blocks_xoz = BC_blocks(1:k,1:k);
            BC_blocks_yoz = BC_blocks(k+1:N,k+1:N);
            
            BC_xoz = obj.build_block_matrix(BC_blocks_xoz);
            BC_yoz = obj.build_block_matrix(BC_blocks_yoz);
            
        end
        
        
        function BC = build_block_matrix(obj, BC_blocks)
            % 构造不规则 block cell 拼接的大矩阵，并自动加上单位阵
            %
            % 输入：
            %   BC_blocks : N×N cell，每个元素是 r_ij × c_ij 的矩阵（可以为空 [])
            %
            % 输出：
            %   BC : 大矩阵 + 单位阵（稀疏格式）
            
            % 找每块的最大尺寸（按行、按列）
            N = size(BC_blocks, 1);
            
            row_sizes = zeros(1, N);
            col_sizes = zeros(1, N);
            for i = 1:N
                for j = 1:N
                    if ~isempty(BC_blocks{i,j})
                        [r,c] = size(BC_blocks{i,j});
                        row_sizes(i) = max(row_sizes(i), r);
                        col_sizes(j) = max(col_sizes(j), c);
                    end
                end
            end
            
            % 先算出最终矩阵大小
            total_rows = sum(row_sizes);
            total_cols = sum(col_sizes);
            
            BigMat = zeros(total_rows, total_cols);
            
            r0 = 0;
            for i = 1:N
                c0 = 0;
                for j = 1:N
                    block = BC_blocks{i,j};
                    if ~isempty(block)
                        [r,c] = size(block);
                        BigMat(r0+1:r0+r, c0+1:c0+c) = block;
                    end
                    c0 = c0 + col_sizes(j);
                end
                r0 = r0 + row_sizes(i);
            end
            
            % 最终一次性拼接成大块矩阵（最快方式）
            % BigMat = cell2mat(BC_blocks);
            n = size(BigMat, 1);     % 或 size(BigMat,2)，因为方阵一样
            sprintf("BC矩阵大小:%dx%d\n",n,n);
            BC = BigMat + sparse(1:n, 1:n, 1, n, n);
        end
        
        
        % 计算出一个BCxx
        function BCxx = splice_BCxx(obj, idx, edge_idx, BC_funcs, BC_loc)
            % 首先拼接每个单独的Q矩阵，Q13_c_h1_h3代表1矩阵与3矩阵相交的c方向分量，横坐标为h1，纵坐标为h3
            idx_case = obj.regions{idx}.impl;       % 当前区域的实例
            edge_case = obj.regions{edge_idx}.impl; % 当前区域邻接的区域的实例
            
            funcss = BC_funcs{idx};  % 这个区域的邻接方程，这里仍然是一个二维矩阵nx6
            
            % 找到非0元素的下标，顺序为c0 c d0 d e f
            row_exists = find(idx_case.coeffs_exists);
            col_exists = find(edge_case.coeffs_exists);
            
            row_has_cd0x = (row_exists(1) == 1);    % 如果有c0一定有d0，一定有这里面的1
            
            rows_len = arrayfun(@(x) (x==1 || x==3).*1 + (x==2 || x==4).*idx_case.H_max + (x>=5).*idx_case.N_max, row_exists);    % 每一行的长度
            cols_len = arrayfun(@(x) (x==1 || x==3).*1 + (x==2 || x==4).*edge_case.H_max + (x>=5).*edge_case.N_max, col_exists);
            
            % 总行数与总列数
            rows_bc = sum(rows_len);
            cols_bc = sum(cols_len);
            
            % 创建一个稀疏矩阵
            BCxx = zeros(rows_bc, cols_bc);
            
            % 当前处理的行，这个函数只会处理一行
            % 这里记录的就是这个邻接区域对应的边界，其中c=B，d=T，e=L，f=R
            edge_bc_loc = BC_loc{idx}(:,edge_idx);  % 对应这个邻接区域的方程的种类，2c 4d 5e 6f
            i = find(row_exists == edge_bc_loc(1));  % 找到这个边界函数对应的下标
            
            % 行数
            if(edge_bc_loc(1) <= 4)  % 是上下边界
                row_hn = idx_case.h;
            else
                row_hn = idx_case.n;
            end
            
            % 这里不是按c d e f的顺序，而是按区域编号顺序
            for j=1:length(col_exists)
                % 找到对应的方程
                func = funcss{edge_bc_loc(2), col_exists(j)};    % 这是这个idx对应的方程
                
                if (col_exists(j) <= 4)
                    col_hn = edge_case.h;
                else
                    col_hn = edge_case.n;
                end
                % 创建一个用于代入数据的矩阵
                [row_idx, col_idx] = ndgrid(1:rows_len(i), 1:cols_len(j));
                
                expr(row_hn, col_hn) = simplifyFraction(rhs(func));
                f = matlabFunction(expr, "Vars", {row_hn, col_hn});
                % Q = arrayfun(@(x,y) (-f(x,y)), row_idx, col_idx);
                
                Q = zeros(size(row_idx));   % 保持矩阵形状
                for ri = 1:size(row_idx,1)
                    for rj = 1:size(row_idx,2)
                        Q(ri,rj) = -f( row_idx(ri,rj), col_idx(ri,rj) );
                    end
                end
                
                % 将这个值送入BCxx中
                % 每处理一个区域，只会有一行的数据，除了c0 d0这样的以外写入2行
                start_row = sum(rows_len(1:i-1)) + 1;
                start_col = sum(cols_len(1:j-1)) + 1;
                end_row = start_row + rows_len(i) - 1;
                end_col = start_col + cols_len(j) - 1;
                BCxx(start_row:end_row, start_col:end_col) = Q;
                
                if edge_bc_loc(1) <= 4 && row_has_cd0x % 如果现在计算的是c或者d，才需要进入这里计算
                    func = funcss{edge_bc_loc(2)-1, col_exists(j)}; % 目前而言，c0/d0方程就是对应c或d的方程的前一个
                    col_idx = 1:cols_len(j);   % 行数为1
                    
                    
                    % row_idx = zeros(1, cols_len(j));
                    % 计算这个矩阵的值
                    % expr = subs(rhs(func), {row_hn,col_hn}, {row_idx, col_idx});
                    % Q1 = -double(expr);
                    expr(col_hn) = simplifyFraction(rhs(func));
                    % pretty(func);
                    f = matlabFunction((expr), "Vars", {col_hn});
                    Q1 = arrayfun(@(x) -f(x), col_idx);
                    
                    % Q1 = zeros(size(col_idx));   % 保持矩阵形状
                    % for ri = 1:size(col_idx,1)
                    %     for rj = 1:size(col_idx,2)
                    %         Q1(ri,rj) = -f(col_idx(ri,rj));
                    %     end
                    % end
                    
                    
                    start_row = start_row - 1;    % 上面的start_row
                    end_row = start_row;    % 这个只需要写一行，写在对应的c或d的上面一行
                    BCxx(start_row:end_row, start_col:end_col) = Q1;
                end
            end
        end
        
        function [ES_xoz, ES_yoz] = splice_ES(obj, ES_funcs)
            % 所有的ES矩阵就是直接拼接，然后转置即可
            N = obj.all_region_num;
            
            ES = [];
            for i=1:N
                if i == obj.region_num_xoz + 1
                    % 如果第一个区域的计算完成了，储存这个区域的
                    ES_xoz = ES';
                    ES = [];
                end
                
                funcs = ES_funcs{i};
                idx_case = obj.regions{i}.impl;       % 当前区域的实例
                
                col_exists = find(idx_case.coeffs_exists);  % 找到当前区域的参数的下标
                
                H_max = idx_case.H_max;
                N_max = idx_case.N_max;
                
                ESx = [];
                
                for j=1:length(col_exists)
                    func_num = col_exists(j);   % 当前读到的方程位置
                    if func_num <= 4   % c或d
                        col_hn = idx_case.h;
                        col_HN_max = H_max;
                    else
                        col_hn = idx_case.n;
                        col_HN_max = N_max;
                    end
                    
                    % 对于每一列，如果real_funcs不存在，那么创建一个0矩阵，否则创建对应矩阵
                    if isempty(funcs{func_num})
                        ESxx = zeros(1,col_HN_max);
                    else
                        func = funcs{func_num};
                        % pretty(func);
                        expr(col_hn) = simplifyFraction(rhs(func));
                        f = matlabFunction(expr, "Vars", {col_hn});
                        
                        if idx_case.ES_regions(i)   % 如果自己这个区域就是有源项
                            % col_idx = 0;
                            ESxx = -double(expr);   % 这里需要加负号
                        else
                            col_idx= 1:col_HN_max;
                            ESxx = arrayfun(@(h) f(h), col_idx);
                            
                            % ESxx = zeros(size(col_idx));   % 保持矩阵形状
                            % for ri = 1:size(col_idx,1)
                            %     for rj = 1:size(col_idx,2)
                            %         ESxx(ri,rj) = f(col_idx(ri,rj));
                            %     end
                            % end
                            
                        end
                        
                    end
                    ESx = [ESx, ESxx];
                    
                end
                ES = [ES, ESx];
            end
            % 剩下的是第二个区域的
            ES_yoz = ES';
            
        end
        
        
        % 用于分割IC矩阵，分割成cell，按c0 c d0 d e f的顺序排列，如果没有则为空
        function [ICs_xoz, ICs_yoz] = split_IC(obj, IC_xoz, IC_yoz)
            all_regions_num = obj.all_region_num;
            len_IC = [obj.len_IC_xoz; obj.len_IC_yoz];
            IC = [IC_xoz; IC_yoz];
            ICs = cell(all_regions_num, 6);
            
            for i=1:all_regions_num
                
                idx_case = obj.regions{i}.impl;     % 当前区域的实例
                % 找到非0元素的下标，顺序为c0 c d0 d e f
                coeffs_exists = idx_case.coeffs_exists;
                H_max = idx_case.H_max;
                N_max = idx_case.N_max;
                
                
                ic_start = sum(len_IC(1:(i-1))) + 1;   % 当前取到的下标的位置
                
                for j=1:length(coeffs_exists)
                    if coeffs_exists(j)
                        ic_end = ic_start + (j==1 || j==3)*1 + (j==2 || j==4)*H_max + (j==5 || j==6)*N_max - 1;
                        ICs{i,j} = IC(ic_start:ic_end);
                        ic_start = ic_end + 1;
                    else
                        ICs{i,j} = [];  % 没有就是空矩阵
                    end
                end
            end
            
            ICs_xoz = ICs(1:obj.region_num_xoz, :);
            ICs_yoz = ICs(obj.region_num_xoz+1:obj.region_num_xoz+obj.region_num_yoz, :);
        end
        
        function cal_IC_lens(obj)
            obj.len_IC_xoz = zeros(obj.region_num_xoz, 1);
            
            for i=1:obj.region_num_xoz
                idx_impl = obj.regions{i}.impl;    % 当前区域的实例
                coeffs_exists = idx_impl.coeffs_exists;
                H_max = idx_impl.H_max;
                N_max = idx_impl.N_max;
                
                obj.len_IC_xoz(i) = coeffs_exists(1) * 1 + coeffs_exists(3) * 1 + ...
                    coeffs_exists(2) .* H_max + coeffs_exists(4) .* H_max + ...
                    coeffs_exists(5) .* H_max + coeffs_exists(6) .* N_max;
            end
            
            obj.len_IC_yoz = zeros(obj.region_num_yoz, 1);
            for i= 1:obj.region_num_yoz
                idx_impl = obj.regions{i+obj.region_num_xoz}.impl;    % 当前区域的实例
                coeffs_exists = idx_impl.coeffs_exists;
                H_max = idx_impl.H_max;
                N_max = idx_impl.N_max;
                
                obj.len_IC_yoz(i) = coeffs_exists(1) * 1 + coeffs_exists(3) * 1 + ...
                    coeffs_exists(2) .* H_max + coeffs_exists(4) .* H_max + ...
                    coeffs_exists(5) .* H_max + coeffs_exists(6) .* N_max;
            end
        end
        
        
        % 计算某个区域某个点集(x0,y0)的Bx与By
        % 需要输入的x0和y0的个数相同，均为行向量（1，n）
        function [Bx, By] = cal_Bx_By_region(obj, plane, ICs, region_num, x0, y0)
            idx_impl = obj.regions{region_num}.impl;    % 当前区域的实例
            if plane == "yoz"
                region_num = region_num - obj.region_num_xoz;   % 给后面的ICs使用
            end
            
            points_num = size(x0, 2);
            
            H_max = idx_impl.H_max;
            N_max = idx_impl.N_max;
            
            Bx_x_expr = idx_impl.B_x_x;    % c d c0 d0
            Bx_y_expr = idx_impl.B_x_y;    % e f
            By_x_expr = idx_impl.B_y_x;    % c d
            By_y_expr = idx_impl.B_y_y;    % e f
            
            c_0x = (idx_impl.c_0x);   % 将 idx_impl.c_0x 转换为符号变量
            c_hx = (idx_impl.c_hx);   % 将 idx_impl.c_hx 转换为符号变量
            d_0x = (idx_impl.d_0x);   % 将 idx_impl.d_0x 转换为符号变量
            d_hx = (idx_impl.d_hx);   % 将 idx_impl.d_hx 转换为符号变量
            e_ny = (idx_impl.e_ny);   % 将 idx_impl.e_ny 转换为符号变量
            f_ny = (idx_impl.f_ny);   % 将 idx_impl.f_ny 转换为符号变量
            
            % 假设 idx_impl 中的成员需要转换为符号变量
            h = sym(idx_impl.h);         % 将 idx_impl.h 转换为符号变量
            n = sym(idx_impl.n);         % 将 idx_impl.n 转换为符号变量
            syms x y real;
            
            hx = repmat((1:H_max)', 1, points_num);
            nx = repmat((1:N_max)', 1, points_num);
            x0h = repmat(x0, H_max, 1);
            y0h = repmat(y0, H_max, 1);
            x0n = repmat(x0, N_max, 1);
            y0n = repmat(y0, N_max, 1);
            
            
            Bx_x_c0 = zeros(1, points_num);
            Bx_x_c = zeros(1, points_num);
            By_x_c = zeros(1, points_num);
            Bx_x_d0 = zeros(1, points_num);
            Bx_x_d = zeros(1, points_num);
            By_x_d = zeros(1, points_num);
            Bx_y_e = zeros(1, points_num);
            By_y_e = zeros(1, points_num);
            Bx_y_f = zeros(1, points_num);
            By_y_f = zeros(1, points_num);
            
            ICn = ICs(region_num,:);    % 取出对应的一行
            for i=1:6
                if ~isempty(ICn{i})
                    switch i
                        case 1
                            Bx_x_c0 = repmat(-ICn{1}, 1, points_num);
                        case 2
                            cx = repmat(ICn{2}, 1, points_num);
                            
                            expr = subs(Bx_x_expr, {c_0x, d_0x, d_hx}, {0, 0, 0});
                            f = matlabFunction(expr, "Vars", {c_hx, h, x, y});
                            % Q = arrayfun(@(c_hx, h, x, y) f(c_hx, h, x, y), cx, hx, x0h, y0h);
                            Q = f(cx, hx, x0h, y0h);
                            Bx_x_c = sum(Q, 1); % 对每一列求和
                            
                            expr = subs(By_x_expr, {c_0x, d_0x, d_hx}, {0, 0, 0});
                            f = matlabFunction(expr, "Vars", {c_hx, h, x, y});
                            % Q = arrayfun(@(c_hx, h, x, y) f(c_hx, h, x, y), cx, hx, x0h, y0h);
                            Q = f(cx, hx, x0h, y0h);
                            By_x_c = sum(Q, 1);
                        case 3
                            Bx_x_d0 = repmat(ICn{3}, 1, points_num);
                        case 4
                            dx = repmat(ICn{4}, 1, points_num);
                            
                            expr = subs(Bx_x_expr, {c_0x, c_hx, d_0x}, {0, 0, 0});
                            f = matlabFunction(expr, "Vars", {d_hx, h, x, y});
                            % Q = arrayfun(@(d_hx, h, x, y) f(d_hx, h, x, y), dx, hx, x0h, y0h);
                            Q = f(dx, hx, x0h, y0h);
                            Bx_x_d = sum(Q, 1);
                            
                            expr = subs(By_x_expr, {c_0x, c_hx, d_0x}, {0, 0, 0});
                            f = matlabFunction(expr, "Vars", {d_hx, h, x, y});
                            % Q = arrayfun(@(d_hx, h, x, y) f(d_hx, h, x, y), dx, hx, x0h, y0h);
                            Q = f(dx, hx, x0h, y0h);
                            By_x_d = sum(Q, 1);
                        case 5
                            ex = repmat(ICn{5}, 1, points_num);
                            
                            expr = subs(Bx_y_expr, f_ny, 0);
                            f = matlabFunction(expr, "Vars", {e_ny, n, x, y});
                            % Q = arrayfun(@(e_ny, n, x, y) f(e_ny, n, x, y), ex, nx, x0n, y0n);
                            Q = f(ex, nx, x0n, y0n);
                            Bx_y_e = sum(Q, 1);
                            
                            expr = subs(By_y_expr, f_ny, 0);
                            f = matlabFunction(expr, "Vars", {e_ny, n, x, y});
                            % Q = arrayfun(@(e_ny, n, x, y) f(e_ny, n, x, y), ex, nx, x0n, y0n);
                            Q = f(ex, nx, x0n, y0n);
                            By_y_e = sum(Q, 1);
                        case 6
                            fx = repmat(ICn{6}, 1, points_num);
                            
                            expr = subs(Bx_y_expr, e_ny, 0);
                            f = matlabFunction(expr, "Vars", {f_ny, n, x, y});
                            % Q = arrayfun(@(f_ny, n, x, y) f(f_ny, n, x, y), fx, nx, x0n, y0n);
                            Q = f(fx, nx, x0n, y0n);
                            Bx_y_f = sum(Q, 1);
                            
                            expr = subs(By_y_expr, e_ny, 0);
                            f = matlabFunction(expr, "Vars", {f_ny, n, x, y});
                            % Q = arrayfun(@(f_ny, n, x, y) f(f_ny, n, x, y), fx, nx, x0n, y0n);
                            Q = f(fx, nx, x0n, y0n);
                            By_y_f = sum(Q, 1);
                    end
                end
            end
            
            
            Bx = Bx_x_c0 + Bx_x_c + Bx_x_d0 + Bx_x_d + Bx_y_e + Bx_y_f;
            By = By_x_c + By_x_d + By_y_e + By_y_f;
            
            if ~isempty(ICn{1}) % 如果存在c0或d0，说明这个区域是一个有源区域，需要加上B_x_P
                B_x_P_v = idx_impl.B_x_P;
                f = matlabFunction(B_x_P_v, "Vars", {y});
                Q = f(y0);
                % Q = arrayfun(@(y) f(y), y0);
                Bx = Bx + Q;
            end
            
            Bx = Bx';
            By = By';
        end
        
        function [Bx, By] = cal_Bx_By(obj, plane, ICs, x0, y0)
            % 注意输入需要是行向量
            if plane == "xoz"
                areas = obj.regions_area(1:obj.region_num_xoz,:);
                [xs, ys, ids] = split_curve_by_rects_by_points(x0, y0, areas);
            elseif plane == "yoz"
                areas = obj.regions_area(obj.region_num_xoz+1:obj.region_num_xoz+obj.region_num_yoz, :);
                [xs, ys, ids] = split_curve_by_rects_by_points(x0, y0, areas);
                ids = ids + obj.region_num_xoz;
            end
            
            segments_lens = length(ids);   % 分区个数
            Bxc = cell(1,segments_lens);
            Byc = cell(1,segments_lens);
            
            parfor i=1:segments_lens
                % 这里只根据ICs计算，和区域编号有关
                [Bxc{i}, Byc{i}] = obj.cal_Bx_By_region(plane, ICs, ids(i) , xs{i}, ys{i});
            end
            
            Bx = vertcat(Bxc{:});
            By = vertcat(Byc{:});
        end
        
        
        function input_current_region(obj, plain, xl, xr, yb, yt, mu_r, J_r)
            % 输入电流区域的坐标与电流大小，线圈匝数
            % 输入的xl,xr,yb,yt分别是左侧x坐标，右侧x坐标，下侧y坐标，上侧y坐标，mu_r指的是这个区域的相对磁导率，一般为1，
            % I_r是电流大小，N_t是线圈匝数
            if plain == "xoz"
                obj.divide_xoz.set_current_regions(xl, xr, yb, yt, mu_r, J_r);
            elseif plain == "yoz"
                obj.divide_yoz.set_current_regions(xl, xr, yb, yt, mu_r, J_r);
            end
        end
        
        function input_calculate_area(obj, xl, xr, yb, yt, mu_r)
            % 输入的xl,xr,yb,yt分别是左侧x坐标，右侧x坐标，下侧y坐标，上侧y坐标，mu_r指的是这个区域的相对磁导率，一般为1
            obj.divide_xoz.set_calculate_area(xl, xr, yb, yt, mu_r);
            obj.divide_yoz.set_calculate_area(xl, xr, yb, yt, mu_r);
        end
        
        
        function pre_process(obj)
            % 输入完所有的区域后，调用这个函数进行预处理
            obj.divide_xoz.divide_regions();
            obj.divide_xoz.findNeighbors();
            obj.divide_xoz.cal_other_info();
            
            % 计算完成后，将数据输入到这个类中
            [regions_area1, obj.region_num_xoz, current_regions1] = obj.divide_xoz.rtn_regions();
            [all_H_max1, all_N_max1] = obj.divide_xoz.rtn_HN_max();
            [all_mu_r1, all_J_r1] = obj.divide_xoz.rtn_mu_J();
            [all_BC_types1, all_casetypes1] = obj.divide_xoz.rtn_types();
            [all_lefts1, all_rights1, all_tops1, all_bottoms1] = obj.divide_xoz.rtn_boundarys_idx(0);
            current_regions_idx1 = obj.divide_xoz.rtn_current_idx(0);
            
            obj.divide_yoz.divide_regions();
            obj.divide_yoz.findNeighbors();
            obj.divide_yoz.cal_other_info();
            
            % 计算完成后，将数据输入到这个类中
            [regions_area2, obj.region_num_yoz, current_regions2] = obj.divide_yoz.rtn_regions();
            [all_H_max2, all_N_max2] = obj.divide_yoz.rtn_HN_max();
            [all_mu_r2, all_J_r2] = obj.divide_yoz.rtn_mu_J();
            [all_BC_types2, all_casetypes2] = obj.divide_yoz.rtn_types();
            [all_lefts2, all_rights2, all_tops2, all_bottoms2] = obj.divide_yoz.rtn_boundarys_idx(obj.region_num_xoz);
            current_regions_idx2 = obj.divide_yoz.rtn_current_idx(obj.region_num_xoz);
            
            obj.regions_area        = [regions_area1; regions_area2];
            obj.current_regions     = [current_regions1; current_regions2];
            obj.all_H_max           = [all_H_max1; all_H_max2];
            obj.all_N_max           = [all_N_max1; all_N_max2];
            obj.all_mu_r            = [all_mu_r1; all_mu_r2];
            obj.all_J_r             = [all_J_r1; all_J_r2];
            obj.all_BC_types        = [all_BC_types1; all_BC_types2];
            obj.all_casetypes       = [all_casetypes1; all_casetypes2];
            obj.all_lefts           = [all_lefts1; all_lefts2];
            obj.all_rights          = [all_rights1; all_rights2];
            obj.all_tops            = [all_tops1; all_tops2];
            obj.all_bottoms         = [all_bottoms1; all_bottoms2];
            obj.current_regions_idx = [current_regions_idx1; current_regions_idx2];
            
            obj.all_region_num = obj.region_num_xoz + obj.region_num_yoz;
        end
        
    end
end