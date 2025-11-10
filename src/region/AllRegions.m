classdef AllRegions < handle
    properties
        regions  % cell array 存储所有区域
        region_num = 7; % 区域数量
        all_H_max = [60; 60; 60; 60; 60; 60; 60];
        all_N_max = [60; 60; 60; 60; 60; 60; 60];
        
        % all_H_max = [300; 300; 107; 107; 43; 21; 21];
        % all_N_max = [60; 60; 268; 268; 268; 268; 268];
        
        % all_H_max = [1; 1; 1; 1; 1; 1; 1];
        % all_N_max = [1; 1; 1; 1; 1; 1; 1];
        
        % all_H_max = [2;2;2;2;2;2;2];
        % all_N_max = [2;2;2;2;2;2;2];
        Jz6 = -1600 .* 5 / (800 .* 1e-6); % Region 6 的电流密度
        Jz7 = 1600 .* 5 / (800 .* 1e-6); % Region 7 的电流密度 (反向)
        
        all_mu_r = [1;1;1;1;1500;1;1];
        all_J_r = [0;0;0;0;0;1e7;-1e7];
        
        len_IC = [60;60;180;180;240;242;242];
    end
    
    methods
        function obj = AllRegions()
            obj.regions = cell(1,7);
        end
        
        function get_all_regions(obj)
            
            obj.set_all_regions();
            disp("区域内方程计算完毕");
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
            ES = obj.splice_ES(ES_funcs);
            disp("计算ES方程完成")
            save("ES.mat", 'ES');
            
            % 开始拼接所有的矩阵
            BC = obj.splice_BC(BC_funcs, BC_loc);
            disp("计算BC方程完成");
            save("BC.mat", 'BC');
            
            writematrix(BC, 'BC_S.xlsx');
            writematrix(ES, 'ES_S.xlsx');
            % IC = BC \ ES;
            IC = lsqr(BC, ES, 1e-6, 1000);
            % IC = pcg(BC, ES, 1e-12, 10000);
            
            save("IC.mat", 'IC');
            
            obj.plot_figures(IC);
            
        end
        
        
        function set_all_regions(obj)
            %% 区域 1
            obj.regions{1} = Region(1, CaseType.BTAir, 0, 0.280, 0, 0.100, ...
                BC_TYPE.BBAA, [3,4,5,6,7],[],[],[],[6,7], obj.all_H_max(1), obj.all_N_max(1), obj.all_mu_r(1), obj.all_J_r(1), obj.region_num);
            obj.regions{1}.get_region_solution_func();
            
            %% 区域 2
            obj.regions{2} = Region(2, CaseType.BTAir, 0, 0.280, 0.140, 0.240, ...
                BC_TYPE.BBAA, [], [3,4,5,6,7],[],[],[6,7], obj.all_H_max(2), obj.all_N_max(2), obj.all_mu_r(2), obj.all_J_r(2), obj.region_num);
            obj.regions{2}.get_region_solution_func();
            
            %% 区域 3
            obj.regions{3} = Region(3, CaseType.NormalAir, 0, 0.100, 0.100, 0.140, ...
                BC_TYPE.AAAA, [2],[1],[],[6],[6,7], obj.all_H_max(3), obj.all_N_max(3), obj.all_mu_r(3), obj.all_J_r(3), obj.region_num);
            obj.regions{3}.get_region_solution_func();
            
            %% 区域 4
            obj.regions{4} = Region(4, CaseType.NormalAir, 0.180, 0.280, 0.100, 0.140, ...
                BC_TYPE.AAAA, [2],[1],[7],[],[6,7], obj.all_H_max(4), obj.all_N_max(4), obj.all_mu_r(4), obj.all_J_r(4), obj.region_num);
            obj.regions{4}.get_region_solution_func();
            
            %% 区域 5
            obj.regions{5} = Region(5, CaseType.NormalAir, 0.120, 0.160, 0.100, 0.140, ...
                BC_TYPE.AAAA, [2],[1],[6],[7],[6,7], obj.all_H_max(5), obj.all_N_max(5), obj.all_mu_r(5), obj.all_J_r(5), obj.region_num);
            obj.regions{5}.get_region_solution_func();
            
            %% 区域6
            obj.regions{6} = Region(6, CaseType.FerriteCurrent, 0.100, 0.120, 0.100, 0.140, ...
                BC_TYPE.AABB, [2],[1],[3],[5],[6,7], obj.all_H_max(6), obj.all_N_max(6), obj.all_mu_r(6), obj.all_J_r(6), obj.region_num);
            obj.regions{6}.get_region_solution_func();
            
            %% 区域7
            obj.regions{7} = Region(7, CaseType.FerriteCurrent, 0.160, 0.180, 0.100, 0.140, ...
                BC_TYPE.AABB, [2],[1],[5],[4],[6,7], obj.all_H_max(7), obj.all_N_max(7), obj.all_mu_r(7), obj.all_J_r(7), obj.region_num);
            obj.regions{7}.get_region_solution_func();
        end
        
        function [BC_funcs, BC_loc, ES_funcs] = cal_all_BCs(obj)
            BC_funcs = cell(1,7);
            BC_loc = cell(1,7);
            ES_funcs = cell(1,7);
            parfor i=1:7    % 并行计算所有边界函数
                [BC_funcs{i}, BC_loc{i}, ES_funcs{i}] = obj.regions{i}.gen_region_coefficient_func(obj.regions);
            end
        end
        
        % 拼接所有的BC
        function BC = splice_BC(obj, BC_funcs, BC_loc)
            
            BC_blocks = cell(obj.region_num, obj.region_num);
            N = obj.region_num;
            sprintf("开始计算BC矩阵，共%d个\n",N);
            parfor i = 1:N
                for j = 1:N
                    if i == j
                        BC_blocks{i,j} = [];    % 暂时置为空
                        continue;
                    end
                    
                    if obj.regions{i}.all_edge_regions(j)   % 有边界相邻
                        BC_blocks{i,j} = obj.splice_BCxx(i, j, BC_funcs, BC_loc);
                    else
                        BC_blocks{i,j} = []; % 或保持空
                    end
                    
                    fprintf("计算完成区域(%d,%d)\n",i,j);
                end
            end
            
            % 找每块的最大尺寸（按行、按列）
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
                % 计算这个矩阵的值
                % 需要先代入数据
                % expr = subs(rhs(func), {row_hn,col_hn}, {row_idx, col_idx});
                % Q = -double(expr);
                expr(row_hn, col_hn) = simplifyFraction(rhs(func));
                f = matlabFunction((expr), "Vars", {row_hn, col_hn});
                Q = arrayfun(@(x,y) -f(x,y), row_idx, col_idx);
                
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
                    % Q1 = -f(row_idx, col_idx);   % 纯数值运算，超级快
                    
                    start_row = start_row - 1;    % 上面的start_row
                    end_row = start_row;    % 这个只需要写一行，写在对应的c或d的上面一行
                    BCxx(start_row:end_row, start_col:end_col) = Q1;
                end
            end
        end
        
        function ES = splice_ES(obj, ES_funcs)
            % 所有的ES矩阵就是直接拼接，然后转置即可
            N = obj.region_num;
            ES = [];
            for i=1:N
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
                            col_idx = 0;
                            ESxx = -double(expr);   % 这里需要加负号
                        else
                            col_idx= 1:col_HN_max;
                            ESxx = arrayfun(@(h) f(h), col_idx);
                        end
                        
                    end
                    ESx = [ESx, ESxx];
                    
                end
                ES = [ES, ESx];
            end
            ES = ES';
        end
        
        % 用于计算总行数或列数
        function out = calc_col_nums(obj, col_nums, h, n)
            % col_nums, h, n 为相同长度向量
            % 输出 out 为同长度数值向量
            
            % 预分配
            out = zeros(size(col_nums));
            
            % ---- case 1：col_nums ≤ 4 且不为1，3 ----
            idx = (col_nums == 2 || col_nums == 4);
            out(idx) = (col_nums(idx)/2).*h(idx) + (col_nums(idx)/2);
            
            % ---- case 2：col_nums ≤ 4 且为1，3 ----
            idx = (col_nums == 1 || col_nums == 3);
            out(idx) = ((col_nums(idx)-1)/2).*h(idx) + 2;
            
            % ---- case 3：col_nums > 4 ----
            idx = (col_nums > 4);
            out(idx) = 2.*h(idx) + 2 + (col_nums(idx)-4).*n(idx);
        end
        
        % 计算某个区域某个点(x0,y0)的Bx与By
        function [Bx, By] = cal_Bx_By(obj, IC, region_num, x0, y0)
            idx_impl = obj.regions{region_num}.impl;    % 当前区域的实例
            
            % 找到非0元素的下标，顺序为c0 c d0 d e f
            coeffs_exists = idx_impl.coeffs_exists;
            
            H_max = idx_impl.H_max;
            N_max = idx_impl.N_max;
            IC_start = sum(obj.len_IC(1:region_num) - obj.len_IC(1)) + 1;
            IC_end = IC_start + obj.len_IC(region_num) - 1;
            ICs = IC(IC_start:IC_end);
            
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
            
            hx = (1:H_max)';
            nx = (1:N_max)';
            x0h = ones(H_max,1) .* x0;
            y0h = ones(H_max,1) .* y0;
            x0n = ones(N_max,1) .* x0;
            y0n = ones(N_max,1) .* y0;
            
            Bx_x_c0 = 0;
            Bx_x_c = 0;
            By_x_c = 0;
            Bx_x_d0 = 0;
            Bx_x_d = 0;
            By_x_d = 0;
            Bx_y_e = 0;
            By_y_e = 0;
            Bx_y_f = 0;
            By_y_f = 0;
            
            
            ICs_idx_s = 1;    % 当前取到的下标
            ICs_idx_e = 1;
            for i=1:length(coeffs_exists)
                if coeffs_exists(i)
                    switch i
                        case 1
                            c0 = ICs(ICs_idx_s);
                            ICs_idx_s = ICs_idx_s + 1;
                            
                            Bx_x_c0 = -c0;
                        case 2
                            ICs_idx_e = ICs_idx_s + H_max -1;
                            cx = ICs(ICs_idx_s:ICs_idx_e);
                            ICs_idx_s = ICs_idx_e;
                            
                            expr = subs(Bx_x_expr, {c_0x, d_0x, d_hx}, {0, 0, 0});
                            f = matlabFunction(expr, "Vars", {c_hx, h, x, y});
                            Q = arrayfun(@(c_hx, h, x, y) f(c_hx, h, x, y), cx, hx, x0h, y0h);
                            Bx_x_c = sum(Q);
                            
                            expr = subs(By_x_expr, {c_0x, d_0x, d_hx}, {0, 0, 0});
                            f = matlabFunction(expr, "Vars", {c_hx, h, x, y});
                            Q = arrayfun(@(c_hx, h, x, y) f(c_hx, h, x, y), cx, hx, x0h, y0h);
                            By_x_c = sum(Q);
                        case 3
                            d0 = ICs(ICs_idx_s);
                            ICs_idx_s = ICs_idx_s + 1;
                            
                            Bx_x_d0 = d0;
                        case 4
                            ICs_idx_e = ICs_idx_s + H_max -1;
                            dx = ICs(ICs_idx_s:ICs_idx_e);
                            ICs_idx_s = ICs_idx_e;
                            
                            expr = subs(Bx_x_expr, {c_0x, c_hx, d_0x}, {0, 0, 0});
                            f = matlabFunction(expr, "Vars", {d_hx, h, x, y});
                            Q = arrayfun(@(d_hx, h, x, y) f(d_hx, h, x, y), dx, hx, x0h, y0h);
                            Bx_x_d = sum(Q);
                            
                            expr = subs(By_x_expr, {c_0x, c_hx, d_0x}, {0, 0, 0});
                            f = matlabFunction(expr, "Vars", {d_hx, h, x, y});
                            Q = arrayfun(@(d_hx, h, x, y) f(d_hx, h, x, y), dx, hx, x0h, y0h);
                            By_x_d = sum(Q);
                        case 5
                            ICs_idx_e = ICs_idx_s + N_max -1;
                            ex = ICs(ICs_idx_s:ICs_idx_e);
                            ICs_idx_s = ICs_idx_e;
                            
                            expr = subs(Bx_y_expr, f_ny, 0);
                            f = matlabFunction(expr, "Vars", {e_ny, n, x, y});
                            Q = arrayfun(@(e_ny, n, x, y) f(e_ny, n, x, y), ex, nx, x0n, y0n);
                            Bx_y_e = sum(Q);
                            
                            expr = subs(By_y_expr, f_ny, 0);
                            f = matlabFunction(expr, "Vars", {e_ny, n, x, y});
                            Q = arrayfun(@(e_ny, n, x, y) f(e_ny, n, x, y), ex, nx, x0n, y0n);
                            By_y_e = sum(Q);
                        case 6
                            ICs_idx_e = ICs_idx_s + N_max -1;
                            fx = ICs(ICs_idx_s:ICs_idx_e);
                            ICs_idx_s = ICs_idx_e;
                            
                            expr = subs(Bx_y_expr, e_ny, 0);
                            f = matlabFunction(expr, "Vars", {f_ny, n, x, y});
                            Q = arrayfun(@(f_ny, n, x, y) f(f_ny, n, x, y), fx, nx, x0n, y0n);
                            Bx_y_f = sum(Q);
                            
                            expr = subs(By_y_expr, e_ny, 0);
                            f = matlabFunction(expr, "Vars", {f_ny, n, x, y});
                            Q = arrayfun(@(f_ny, n, x, y) f(f_ny, n, x, y), fx, nx, x0n, y0n);
                            By_y_f = sum(Q);
                    end
                end
            end
            
            
            
            Bx = Bx_x_c0 + Bx_x_c + Bx_x_d0 + Bx_x_d + Bx_y_e + Bx_y_f;
            By = By_x_c + By_x_d + By_y_e + By_y_f;
            
            if coeffs_exists(1) % 如果存在c0或d0，说明这个区域是一个有源区域，需要加上B_x_P
                B_x_P_v = subs(idx_impl.B_x_P, y, y0);
                Bx = Bx + double(B_x_P_v);
            end
        end
        
        function plot_figures(obj, IC)
            points = 200;
            x0 = 0.14 .* ones(1,points);
            y0 = linspace(0,0.240,points);
            Bx = zeros(1,points);
            By = zeros(1,points);
            
            for i=1:points
                if y0(i) > 0 && y0(i) <= 0.1
                    [Bx(i), By(i)] = obj.cal_Bx_By(IC, 1, x0(i), y0(i));
                elseif y0(i) > 0.1 && y0(i) <= 0.14
                    [Bx(i), By(i)] = obj.cal_Bx_By(IC, 5, x0(i), y0(i));
                elseif y0(i) > 0.14 && y0(i) <= 0.24
                    [Bx(i), By(i)] = obj.cal_Bx_By(IC, 2, x0(i), y0(i));
                end
            end
            
            figure;
            plot(y0,Bx);    hold on;
            plot(y0,By);
            grid on;
        end
        
    end
end