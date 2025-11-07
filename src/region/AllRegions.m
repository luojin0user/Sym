classdef AllRegions < handle
    properties
        regions  % cell array 存储所有区域
        region_num = 7; % 区域数量
        all_H_max = [60; 60; 60; 60; 60; 60; 60];
        all_N_max = [60; 60; 60; 60; 60; 60; 60];
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
            [BC_funcs, BC_loc] = obj.cal_all_BCs();
            disp("计算系数方程完成");
            % 注意这里BCs的顺序
            % 每个BCs{i}代表对应区域i的所有边界方程
            % BCs{i}{j,k}代表当前区域按顺序的{j,k}个方程，例如{2,3}表示区域的第二个临界区域的第三个方程（含有e的方程）
            % 区域内的方程的每一行代表一个邻接区域，按“上下左右”的顺序排列，如上边界有多个就是“上上上上下左右”
            % 每一行的每个元素代表这个邻接区域的方程的c0 c d0 d e f分量，如果为{}则没有该分量，直接跳过
            
            % 开始拼接所有的矩阵
            BC = obj.splice_BC(BC_funcs, BC_loc);
            disp("计算BC方程完成");
            % disp(BC);
            
            save("BC.mat", 'BC');
            
            %{
            %% 符号代入计算
            eq_main = BC1{4, 1};
            % eq_main_means = M3{4}{1};
            
            % 将 A_z_3 代入 d_hx_3_eq
            % eq_final = subs(eq_main, lhs(eq_sub1), rhs(eq_sub1));
            eq_final = subs(eq_main, obj.regions{2}.impl.h, 1);
            eq_final = subs(eq_final, obj.regions{2}.impl.n, 1);
            eq_final = subs(eq_final, obj.regions{6}.impl.n, 1);
            eq_final = subs(eq_final, obj.regions{6}.impl.h, 1);
            % 去除 Hold，自动求积分
            eq_final = release(rhs(eq_final));
            
            eq_final_value = -double(eq_final); % 系数都需要取负号
            
            disp('代入后的符号方程 eqc:');
            pretty(eq_final);
            disp('代入后的符号方程 eqc 数值:');
            disp(eq_final_value);
            %disp("对应的系数：");
            %disp(eq_main_means);
            %}
        end
        
        
        function set_all_regions(obj)
            %% 区域 1
            obj.regions{1} = Region(1, CaseType.BTAir, 0, 280, 0, 100, ...
                BC_TYPE.BBAA, [3,4,5,6,7],[],[],[], obj.all_H_max(1), obj.all_N_max(1), obj.region_num);
            K1 = obj.regions{1}.get_region_solution_func();
            %pretty(K1{1});  % MATLAB 版本的 pprint
            
            %% 区域 2
            obj.regions{2} = Region(2, CaseType.BTAir, 0, 280, 140, 240, ...
                BC_TYPE.BBAA, [], [3,4,5,6,7],[],[], obj.all_H_max(2), obj.all_N_max(2), obj.region_num);
            K2 = obj.regions{2}.get_region_solution_func();
            
            %pretty(K2{1});  % MATLAB 版本的 pprint
            
            %% 区域 3
            obj.regions{3} = Region(3, CaseType.NormalAir, 0, 100, 100, 140, ...
                BC_TYPE.AAAA, [2],[1],[],[6], obj.all_H_max(3), obj.all_N_max(3), obj.region_num);
            K3 = obj.regions{3}.get_region_solution_func();
            
            %pretty(K3{1});  % MATLAB 版本的 pprint
            
            %% 区域 4
            obj.regions{4} = Region(4, CaseType.NormalAir, 180, 280, 100, 140, ...
                BC_TYPE.AAAA, [2],[1],[7],[], obj.all_H_max(4), obj.all_N_max(4), obj.region_num);
            K4 = obj.regions{4}.get_region_solution_func();
            
            %pretty(K4{1});  % MATLAB 版本的 pprint
            
            %% 区域 5
            obj.regions{5} = Region(3, CaseType.NormalAir, 120, 160, 100, 140, ...
                BC_TYPE.AAAA, [2],[1],[6],[7], obj.all_H_max(5), obj.all_N_max(5), obj.region_num);
            K5 = obj.regions{5}.get_region_solution_func();
            
            %pretty(K5{1});  % MATLAB 版本的 pprint
            
            %% 区域6
            obj.regions{6} = Region(6, CaseType.FerriteCurrent, 100, 120, 100, 140, ...
                BC_TYPE.AABB, [2],[1],[3],[5], obj.all_H_max(6), obj.all_N_max(6), obj.region_num);
            K6 = obj.regions{6}.get_region_solution_func();
            % pretty(K6{1});  % MATLAB 版本的 pprint
            
            %% 区域7
            obj.regions{7} = Region(7, CaseType.FerriteCurrent, 160, 180, 100, 140, ...
                BC_TYPE.AABB, [2],[1],[5],[4], obj.all_H_max(7), obj.all_N_max(7), obj.region_num);
            K6 = obj.regions{7}.get_region_solution_func();
            %pretty(K6{1});  % MATLAB 版本的 pprint
        end
        
        function [BC_funcs, BC_loc] = cal_all_BCs(obj)
            BC_funcs = cell(1,7);
            BC_loc = cell(1,7);
            regions_tmp = obj.regions;
            parfor i=1:7    % 并行计算所有边界函数
                [BC_funcs{i}, BC_loc{i}] = regions_tmp{i}.gen_region_coefficient_func(regions_tmp);
            end
        end
        
        % 拼接所有的BC
        function BC = splice_BC(obj, BC_funcs, BC_loc)
            
            BC_blocks = cell(obj.region_num, obj.region_num);
            N = obj.region_num;
            fprintf("开始计算BC矩阵，共%d个\n",N);
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
                expr(row_hn, col_hn) = release(rhs(func));
                f = matlabFunction(expr, "Vars", {row_hn, col_hn});
                Q = -f(row_idx, col_idx);   % 纯数值运算，超级快
                
                % 将这个值送入BCxx中
                % 每处理一个区域，只会有一行的数据，除了c0 d0这样的以外写入2行
                start_row = sum(rows_len(1:i-1)) + 1;
                start_col = sum(cols_len(1:j-1)) + 1;
                end_row = start_row + rows_len(i) - 1;
                end_col = start_col + cols_len(j) - 1;
                BCxx(start_row:end_row, start_col:end_col) = Q;
                
                if edge_bc_loc(1) <= 4 && row_has_cd0x % 如果现在计算的是c或者d，才需要进入这里计算
                    func = funcss{edge_bc_loc(2)-1, col_exists(j)}; % 目前而言，c0/d0方程就是对应c或d的方程的前一个
                    [~, col_idx] = ndgrid(1, 1:cols_len(j));   % 行数为1
                    row_idx = zeros(1, cols_len(j));
                    % 计算这个矩阵的值
                    % expr = subs(rhs(func), {row_hn,col_hn}, {row_idx, col_idx});
                    % Q1 = -double(expr);
                    expr(row_hn, col_hn) = release(rhs(func));
                    f = matlabFunction(expr, "Vars", {row_hn, col_hn});
                    Q1 = -f(row_idx, col_idx);   % 纯数值运算，超级快
                    
                    end_row = start_row - 1;    % 上面的start_row
                    start_row = end_row;    % 这个只需要写一行，写在对应的c或d的上面一行
                    BCxx(start_row:end_row, start_col:end_col) = Q1;
                end
            end
        end
    end
end