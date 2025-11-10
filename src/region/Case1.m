classdef Case1 < BasicCase
    methods
        function obj = Case1(idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max, mu_r)
            % 调用父类构造函数
            obj@BasicCase(idx, xl, xr, yl, yt, H_max, N_max, mu_r);
            obj.num_coeffs = 4;
            obj.apply_boundaries(Ln, Rn, Tn, Bn);
            
            obj.coeffs_exists = [0; ~Bn; 0; ~Tn; ~Ln; ~Rn];
        end
        
        function gen_solution_func(obj)
            syms x y
            % 总的 A_z
            A_z(x,y) = obj.A_zx_expr + obj.A_zy_expr;
            
            % 计算磁场分量
            obj.B_x_x = diff(obj.A_zx_expr , y);
            obj.B_x_y = diff(obj.A_zy_expr , y);
            obj.B_y_x = -diff(obj.A_zx_expr , x);
            obj.B_y_y = -diff(obj.A_zy_expr , x);
            
            B_xx = obj.B_x_x + obj.B_x_y;
            B_xy = obj.B_y_x + obj.B_y_y;
            
            % 构造符号方程
            suffix = ['_' num2str(obj.idx)];
            obj.eq_A_z = symfun(sym(['A_z' suffix]), [x y]) == A_z;
            obj.eq_B_x = symfun(sym(['B_x' suffix]), [x y]) == (B_xx + B_xy);
        end
        
        % 设置系数方程
        function gen_coefficient_func(obj)
            syms x y
            suffix = ['_' num2str(obj.idx)];
            % 首先需要知道传入的边界方程，然后直接带入这个方程，这样子就是直接将边界方程带入下面的等式中
            % 也就是经过这个计算之后就直接确定了最后的参数方程
            % 计算 c_hx, d_hx, e_ny, f_ny 的积分表达式
            % 这里需要注意，传入的方程 e.g. obj.B_funcs 是一个cell数组，其中的每一项是这个边界的另一边的函数的c d e f系数，对应 obj.B_coeffs ，所以需要对这个边界的4个系数分别计算，也就最多会产生16个系数方程
            obj.eq_c_hx = cell(1,6);
            obj.eq_d_hx = cell(1,6);
            obj.eq_e_ny = cell(1,6);
            obj.eq_f_ny = cell(1,6);
            
            obj.eq_ES = cell(1,6);
            
            es_c_expr = 0;
            es_d_expr = 0;
            es_e_expr = 0;
            es_f_expr = 0;
            
            
            func_num = 1;
            
            % c0 d0为空
            obj.eq_c0x = cell(1,6);
            func_num = func_num + 1;
            rowk = 0;
            rowES = 1;
            if ~obj.Bn
                for i = 1:length(obj.B_funcs)
                    row = ceil(i/6);
                    col = mod(i+5,6)+1;
                    if isempty(obj.B_funcs{i})
                        obj.eq_c_hx{row, col} = [];
                        continue; % 跳过为零的函数
                    end
                    % 对于分段函数的积分上下限，是对应的邻接区域的上下限，而不是当前区域的上下限
                    % 首先找到对应临界区域
                    bottom_idx = obj.bottoms(row);
                    bottom_i = obj.all_regions{bottom_idx}.impl;  % 对应的top_i的对象的实现
                    % 注意这里的积分限，是自己边和邻接边的相交位置
                    int_start = max(bottom_i.xl, obj.xl);
                    int_end = min(bottom_i.xr, obj.xr);
                    % [int_start, int_end] = obj.find_intersection(bottom_i.xl, bottom_i.xr, obj.xl, obj.xr);
                    c_hx_expr(x,y) = (2 / obj.tau_x) * int(obj.B_funcs{i} * sin(obj.beta_h * (x - obj.xl)), x, int_start, int_end, Hold=true);
                    obj.eq_c_hx{row, col} = symfun(sym(['c_hx' suffix]), [x,y]) == (c_hx_expr);
                    
                    if row ~= rowk
                        obj.BCfuncs_loc_map(:,bottom_idx) = [2,func_num];
                        func_num = func_num + 1;
                        rowk = row;
                        
                        % 如果当前这个临界区域是有源区域
                        if(obj.ES_regions(bottom_idx) == true)
                            % 找到对应的方程，可能是求导的，可能是原方程
                            es_expr = (2 / obj.tau_x) * int(obj.B_ESfuncs(rowES) * sin(obj.beta_h * (x - obj.xl)), x, int_start, int_end, Hold=true);
                            es_c_expr = es_c_expr + es_expr;    % 在一个边界的所有ES相加
                            rowES = rowES + 1;
                        end
                    end
                end
            else
                func_num = func_num + 1;
            end
            
            obj.eq_d0x = cell(1,6);
            func_num = func_num + 1;
            
            rowk = 0;
            rowES = 1;
            if ~obj.Tn
                for i = 1:length(obj.T_funcs)
                    row = ceil(i/6);
                    col = mod(i+5,6)+1;
                    if isempty(obj.T_funcs{i})
                        obj.eq_d_hx{row, col} = [];
                        continue; % 跳过为零的函数
                    end
                    top_idx = obj.tops(row);
                    top_i = obj.all_regions{top_idx}.impl;  % 对应的top_i的对象的实现
                    int_start = max(top_i.xl, obj.xl);
                    int_end = min(top_i.xr, obj.xr);
                    % [int_start, int_end] = obj.find_intersection(top_i.xl, top_i.xr, obj.xl, obj.xr);
                    d_hx_expr(x,y) = (2 / obj.tau_x) * int(obj.T_funcs{i} * sin(obj.beta_h * (x - obj.xl)), x, int_start, int_end, Hold=true);
                    obj.eq_d_hx{row, col} = symfun(sym(['d_hx' suffix]), [x,y]) ==(d_hx_expr);
                    
                    if row ~= rowk
                        obj.BCfuncs_loc_map(:,top_idx) = [4,func_num];
                        func_num = func_num + 1;
                        rowk = row;
                        
                        % 如果当前这个临界区域是有源区域
                        if obj.ES_regions(top_idx)
                            % 找到对应的方程，可能是求导的，可能是原方程
                            es_expr = (2 / obj.tau_x) * int(obj.T_ESfuncs(rowES) * sin(obj.beta_h * (x - obj.xl)), x, int_start, int_end, Hold=true);
                            es_d_expr = es_d_expr + es_expr;    % 在一个边界的所有ES相加
                            rowES = rowES + 1;
                        end
                    end
                end
            else
                func_num = func_num + 1;
            end
            
            if ~obj.Ln
                for i = 1:length(obj.L_funcs)
                    if isempty(obj.L_funcs{i})
                        obj.eq_e_ny{i} = [];
                        continue; % 跳过为零的函数
                    end
                    e_ny_expr(x,y) = (2 / obj.tau_y) * int(obj.L_funcs{i} * sin(obj.lambda_n * (y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y, Hold=true);
                    obj.eq_e_ny{i} = symfun(sym(['e_ny' suffix]), [x,y]) == (e_ny_expr);
                end
                
                left_idx = obj.lefts(1);
                % 如果当前这个临界区域是有源区域
                if obj.ES_regions(left_idx)
                    % 找到对应的方程，可能是求导的，可能是原方程
                    es_expr = (2 / obj.tau_y) * int(obj.L_ESfuncs(1) * sin(obj.lambda_n * (y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y, Hold=true);
                    es_e_expr = es_e_expr + es_expr;    % 在一个边界的所有ES相加
                end
                
                obj.BCfuncs_loc_map(:,left_idx) = [5,func_num];
                func_num = func_num + 1;
            else
                func_num = func_num + 1;
            end
            
            if ~obj.Rn
                for i = 1:length(obj.R_funcs)
                    if isempty(obj.R_funcs{i})
                        obj.eq_f_ny{i} = [];
                        continue; % 跳过为零的函数
                    end
                    f_ny_expr(x,y) = (2 / obj.tau_y) * int(obj.R_funcs{i} * sin(obj.lambda_n * (y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y, Hold=true);
                    obj.eq_f_ny{i} = symfun(sym(['f_ny' suffix]), [x,y]) == (f_ny_expr);
                end
                
                right_idx = obj.rights(1);
                % 如果当前这个临界区域是有源区域
                if obj.ES_regions(right_idx)
                    % 找到对应的方程，可能是求导的，可能是原方程
                    es_expr = (2 / obj.tau_y) * int(obj.R_ESfuncs(1) * sin(obj.lambda_n * (y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y, Hold=true);
                    es_f_expr = es_f_expr + es_expr;    % 在一个边界的所有ES相加
                end
                
                obj.BCfuncs_loc_map(:,right_idx) = [6,func_num];
                func_num = func_num + 1;
            else
                func_num = func_num + 1;
            end
            
            obj.eq_ES{1} = [];
            obj.eq_ES{3} = [];
            
            % c_ES
            if es_c_expr ~= 0
                obj.eq_ES{2} = symfun(sym(['c_ES' suffix]), [x, y]) == (es_c_expr);
            else
                obj.eq_ES{2} = [];   % 空表示不用创建
            end
            
            % d_ES
            if es_d_expr ~= 0
                obj.eq_ES{4} = symfun(sym(['d_ES' suffix]), [x, y]) == (es_d_expr);
            else
                obj.eq_ES{4} = [];
            end
            
            % e_ES
            if es_e_expr ~= 0
                obj.eq_ES{5} = symfun(sym(['e_ES' suffix]), [x, y]) == (es_e_expr);
            else
                obj.eq_ES{5} = [];
            end
            
            % f_ES
            if es_f_expr ~= 0
                obj.eq_ES{6} = symfun(sym(['f_ES' suffix]), [x, y]) == (es_f_expr);
            else
                obj.eq_ES{6} = [];
            end
            
            
            % 这里直接改变了这个区域的边界方程情况，无需再进行返回值传递
            
            % regions = {obj.eq_c_hx, obj.eq_d_hx, obj.eq_e_ny, obj.eq_f_ny};
        end
    end
end
