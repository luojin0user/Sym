classdef Case2 < BasicCase
    properties
        J_z
        A_y_P
        B_x_P
        B_y_P
    end
    
    methods
        function obj = Case2(idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max, mu_r, J_r)
            % 调用父类构造函数
            obj@BasicCase(idx, xl, xr, yl, yt, H_max, N_max, mu_r);
            % 应用边界条件
            obj.num_coeffs = 6;
            obj.apply_boundaries(Ln, Rn, Tn, Bn);
            
            suffix = ['_' num2str(idx)];
            
            % obj.J_z = sym(['J_z' suffix], 'real');
            obj.J_z = J_r;
            obj.c_0x = sym(['c_0x' suffix], 'real');
            obj.d_0x = sym(['d_0x' suffix], 'real');
            
            % 线性项
            syms x y real
            obj.A_y_P = 0.5 * obj.mu_0 * obj.mu_r * obj.J_z * y^2;
            obj.B_x_P = diff(obj.A_y_P, y);
            obj.B_y_P = -diff(obj.A_y_P, x);
            
            obj.eq_A_y_P = symfun(sym(['A_y_P' suffix]), [x y]) == obj.A_y_P;
            obj.eq_B_x_P = symfun(sym(['B_x_P' suffix]), [x y]) == obj.B_x_P;
            obj.eq_B_y_P = symfun(sym(['B_y_P' suffix]), [x y]) == obj.B_y_P;
            
            obj.coeffs_exists = [~Bn; ~Bn; ~Tn; ~Tn; ~Ln; ~Rn];
        end
        
        function regions = gen_solution_func(obj)
            syms x y
            
            A_z = obj.A_zx_expr + obj.A_zy_expr;
            
            % 磁场分量
            obj.B_x_x = diff(obj.A_zx_expr , y);
            obj.B_x_y = diff(obj.A_zy_expr , y);
            obj.B_y_x = -diff(obj.A_zx_expr , x);
            obj.B_y_y = -diff(obj.A_zy_expr , x);
            
            B_xx = obj.B_x_x + obj.B_x_y;
            B_xy = obj.B_y_x + obj.B_y_y;
            
            % 符号方程
            suffix = ['_' num2str(obj.idx)];
            obj.eq_A_z = symfun(sym(['A_z' suffix]), [x y]) ==  A_z;
            obj.eq_B_x = symfun(sym(['B_x' suffix]), [x y]) == (B_xx + B_xy);
            
            % 返回 cell array
            % regions = {obj.impl.eq_A_z, obj.impl.eq_B_x};
        end
        
        % 设置系数方程
        function gen_coefficient_func(obj)
            syms x y
            suffix = ['_' num2str(obj.idx)];
            
            obj.eq_c0x = cell(1,6);
            obj.eq_d0x = cell(1,6);
            obj.eq_c_hx = cell(1,6);
            obj.eq_d_hx = cell(1,6);
            obj.eq_e_ny = cell(1,6);
            obj.eq_f_ny = cell(1,6);
            
            obj.eq_ES = cell(1,6);
            
            func_num = 1;
            
            if ~obj.Bn
                for i = 1:length(obj.B_funcs)
                    if isempty(obj.B_funcs{i})
                        obj.eq_c_hx{i} = [];
                        obj.eq_c0x{i} = [];
                        continue; % 跳过为零的函数
                    end
                    c_hx_expr(x,y) = (2/obj.tau_x) * int(obj.B_funcs{i} .* cos(obj.beta_h*(x - obj.xl)), x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
                    obj.eq_c_hx{i} = symfun(sym(['c_hx' suffix]), [x,y]) == c_hx_expr;
                    
                    c_0x_expr(x,y) = (1/obj.tau_x) * int((1/obj.tau_y) * obj.B_funcs{i}, x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
                    obj.eq_c0x{i} = symfun(sym(['c_0x' suffix]), [x,y]) == c_0x_expr;
                    
                end
                
                bottom_idx = obj.bottoms(1);
                % 如果当前这个临界区域是有源区域
                if obj.ES_regions(bottom_idx)
                    % 找到对应的方程，可能是求导的，可能是原方程
                    es_expr = (1/obj.tau_x) * int((1/obj.tau_y) * obj.B_ESfuncs(1), x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
                    es_c_expr = es_c_expr + es_expr;    % 在一个边界的所有ES相加
                end
                
                func_num = func_num + 1;
                obj.BCfuncs_loc_map(:,bottom_idx) = [2,func_num];
                func_num = func_num + 1;
            else
                func_num = func_num + 1;
            end
            
            if ~obj.Tn
                for i = 1:length(obj.T_funcs)
                    if isempty(obj.T_funcs{i})
                        obj.eq_d_hx{i} = [];
                        obj.eq_d0x{i} = [];
                        continue; % 跳过为零的函数
                    end
                    d_hx_expr(x,y) = (2/obj.tau_x) * int(obj.T_funcs{i} .* cos(obj.beta_h*(x - obj.xl)), x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
                    obj.eq_d_hx{i} = symfun(sym(['d_hx' suffix]), [x,y]) == d_hx_expr;
                    
                    d_0x_expr(x,y) = (1/obj.tau_x) * int((1/obj.tau_y) * obj.T_funcs{i}, x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
                    obj.eq_d0x{i} = symfun(sym(['d_0x' suffix]), [x,y]) == d_0x_expr;
                    
                end
                
                top_idx = obj.tops(1);
                % 如果当前这个临界区域是有源区域
                if obj.ES_regions(top_idx)
                    % 找到对应的方程，可能是求导的，可能是原方程
                    es_expr = (1/obj.tau_x) * int((1/obj.tau_y) * obj.T_ESfuncs(1), x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
                    es_d_expr = es_d_expr + es_expr;    % 在一个边界的所有ES相加
                end
                
                func_num = func_num + 1;
                obj.BCfuncs_loc_map(:,top_idx) = [4,func_num];
                func_num = func_num + 1;
            else
                func_num = func_num + 1;
            end
            
            if ~obj.Ln
                for i = 1:length(obj.L_funcs)
                    if isempty(obj.L_funcs{i})
                        obj.eq_e_ny{i} = [];
                        continue; % 跳过为零的函数
                    end
                    e_ny_expr(x,y) = (2/obj.tau_y) * int(obj.L_funcs{i} .* sin(obj.lambda_n*(y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y,'Hold',true);
                    obj.eq_e_ny{i} = symfun(sym(['e_ny' suffix]), [x,y]) == e_ny_expr;
                end
                
                left_idx = obj.lefts(1);
                % 如果当前这个临界区域是有源区域
                if obj.ES_regions(left_idx)
                    % 找到对应的方程，可能是求导的，可能是原方程
                    es_expr = (2/obj.tau_y) * int(obj.L_ESfuncs(1) .* sin(obj.lambda_n*(y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y,'Hold',true);
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
                    f_ny_expr(x,y) = (2/obj.tau_y) * int(obj.R_funcs{i} .* sin(obj.lambda_n*(y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y,'Hold',true);
                    obj.eq_f_ny{i} = symfun(sym(['f_ny' suffix]), [x,y]) == f_ny_expr;
                end
                
                right_idx = obj.rights(1);
                % 如果当前这个临界区域是有源区域
                if obj.ES_regions(right_idx)
                    % 找到对应的方程，可能是求导的，可能是原方程
                    es_expr = (2/obj.tau_y) * int(obj.R_ESfuncs(1) .* sin(obj.lambda_n*(y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y,'Hold',true);
                    es_f_expr = es_f_expr + es_expr;    % 在一个边界的所有ES相加
                end
                
                obj.BCfuncs_loc_map(:,right_idx) = [6,func_num];
                func_num = func_num + 1;
            else
                func_num = func_num + 1;
            end
            
            % c_ES c0
            es_c_expr(x,y) = (1/obj.tau_x) * int((1/obj.tau_y) * obj.A_y_P, x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
            es_c_expr = subs(es_c_expr, y, obj.yl);
            obj.eq_ES{1} = symfun(sym(['c_ES' suffix]), [x, y]) == es_c_expr;
            
            % d_ES
            es_d_expr(x,y) = (1/obj.tau_x) * int((1/obj.tau_y) * obj.A_y_P, x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
            es_d_expr = subs(es_d_expr, y, obj.yt);
            obj.eq_ES{3} = symfun(sym(['d_ES' suffix]), [x, y]) == es_d_expr;
            
            
            obj.eq_ES{2} = [];
            obj.eq_ES{4} = [];
            obj.eq_ES{5} = [];
            obj.eq_ES{6} = [];
            
        end
    end
end
