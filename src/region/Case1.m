classdef Case1 < BasicCase
    methods
        function obj = Case1(idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max)
            % 调用父类构造函数
            obj@BasicCase(idx, xl, xr, yl, yt, H_max, N_max);
            obj.apply_boundaries(Ln, Rn, Tn, Bn);
        end
        
        function gen_solution_func(obj)
            syms x y
            % 应用边界条件
            
            % 构造符号求和
            % 沿 x 方向的 A_z 分量
            
            %{
            A_zx_expr(x,y) = symsum( ...
                ( (obj.c_hx / obj.beta_h) * sinh(obj.beta_h * (obj.yt - y)) / sinh(obj.beta_h * obj.tau_y) ...
                + (obj.d_hx / obj.beta_h) * sinh(obj.beta_h * (y - obj.yl)) / sinh(obj.beta_h * obj.tau_y) ) ...
                * sin(obj.beta_h * (x - obj.xl)), ...
                obj.h, 1, obj.H_max);
            
            % 沿 y 方向的 A_z 分量
            A_zy_expr(x,y) = symsum( ...
                ( (obj.e_ny / obj.lambda_n) * sinh(obj.lambda_n * (obj.xr - x)) / sinh(obj.lambda_n * obj.tau_x) ...
                + (obj.f_ny / obj.lambda_n) * sinh(obj.lambda_n * (x - obj.xl)) / sinh(obj.lambda_n * obj.tau_x) ) ...
                * sin(obj.lambda_n * (y - obj.yl)), ...
                obj.n, 1, obj.N_max);
            %}
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
            
            for i = 1:length(obj.B_funcs)
                if obj.B_funcs{i} == 0
                    obj.eq_c_hx{row, ceil(i/row)} = 0;
                    continue; % 跳过为零的函数
                end
                c_hx_expr(x,y) = (2 / obj.tau_x) * int(obj.B_funcs{i} * sin(obj.beta_h * (x - obj.xl)), x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
                
                obj.eq_c_hx{ceil(i/6), mod(i+5,6)} = symfun(sym(['c_hx' suffix]), [x,y]) == c_hx_expr;
            end
            
            for i = 1:length(obj.T_funcs)
                if obj.T_funcs{i} == 0
                    obj.eq_d_hx{row, ceil(i/row)} = 0;
                    continue; % 跳过为零的函数
                end
                d_hx_expr(x,y) = (2 / obj.tau_x) * int(obj.T_funcs{i} * sin(obj.beta_h * (x - obj.xl)), x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
                obj.eq_d_hx{ceil(i/6), mod(i+5,6)+1} = symfun(sym(['d_hx' suffix]), [x,y]) == d_hx_expr;
            end
            
            for i = 1:length(obj.L_funcs)
                if obj.L_funcs{i} == 0
                    obj.eq_e_ny{row, ceil(i/row)} = 0;
                    continue; % 跳过为零的函数
                end
                e_ny_expr(x,y) = (2 / obj.tau_y) * int(obj.L_funcs{i} * sin(obj.lambda_n * (y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y,'Hold',true);
                obj.eq_e_ny{i} = symfun(sym(['e_ny' suffix]), [x,y]) == e_ny_expr;
            end
            
            for i = 1:length(obj.R_funcs)
                if obj.L_funcs{i} == 0
                    obj.eq_f_ny{row, ceil(i/row)} = 0;
                    continue; % 跳过为零的函数
                end
                f_ny_expr(x,y) = (2 / obj.tau_y) * int(obj.R_funcs{i} * sin(obj.lambda_n * (y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y,'Hold',true);
                obj.eq_f_ny{i} = symfun(sym(['f_ny' suffix]), [x,y]) == f_ny_expr;
            end
            
            % 这里直接改变了这个区域的边界方程情况，无需再进行返回值传递
            
            % regions = {obj.eq_c_hx, obj.eq_d_hx, obj.eq_e_ny, obj.eq_f_ny};
        end
    end
end
