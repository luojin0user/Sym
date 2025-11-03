classdef Case2 < BasicCase
    properties
        J_z
        A_y_P
    end
    
    methods
        function obj = Case2(idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max)
            % 调用父类构造函数
            obj@BasicCase(idx, xl, xr, yl, yt, H_max, N_max);
            % 应用边界条件
            obj.apply_boundaries(Ln, Rn, Tn, Bn);
            
            suffix = ['_' num2str(idx)];
            syms(['J_z' suffix], ['c_0x' suffix], ['d_0x' suffix], 'real')
            
            obj.J_z = eval(['J_z' suffix]);
            obj.c_0x = eval(['c_0x' suffix]);
            obj.d_0x = eval(['d_0x' suffix]);
            
            % 线性项
            syms y real
            obj.A_y_P = 0.5 * obj.mu_0 * obj.mu_r * obj.J_z * y^2;
        end
        
        function regions = gen_solution_func(obj)
            syms x y
            
            
            % 构造符号求和
            %{
            % 沿 x 方向 A_z
            A_zx_expr(x,y) = (obj.yt - y) * obj.c_0x + (y - obj.yl) * obj.d_0x + ...
                symsum( ...
                ((obj.c_hx / obj.beta_h) * sinh(obj.beta_h * (obj.yt - y)) / sinh(obj.beta_h * obj.tau_y) + ...
                (obj.d_hx / obj.beta_h) * sinh(obj.beta_h * (y - obj.yl)) / sinh(obj.beta_h * obj.tau_y)) .* ...
                cos(obj.beta_h * (x - obj.xl)), ...
                obj.h, 1, obj.H_max);
            
            % 沿 y 方向 A_z
            A_zy_expr(x,y) = -symsum( ...
                ((obj.e_ny / obj.lambda_n) .* cosh(obj.lambda_n * (x - obj.xl)) ./ sinh(obj.lambda_n * obj.tau_x) - ...
                (obj.f_ny / obj.lambda_n) .* cosh(obj.lambda_n * (obj.xr - x)) ./ sinh(obj.lambda_n * obj.tau_x)) .* ...
                sin(obj.lambda_n * (y - obj.yl)), ...
                obj.n, 1, obj.N_max);
            %}
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
            obj.eq_c0x = cell(1,4);
            obj.eq_d0x = cell(1,4);
            obj.eq_c_hx = cell(1,4);
            obj.eq_d_hx = cell(1,4);
            obj.eq_e_ny = cell(1,4);
            obj.eq_f_ny = cell(1,4);
            
            % 符号方程封装
            suffix = ['_' num2str(obj.idx)];
            
            for i = 1:length(obj.B_funcs)
                c_hx_expr(x,y) = (2/obj.tau_x) * int(obj.B_funcs{i} .* cos(obj.beta_h*(x - obj.xl)), x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
                obj.eq_c_hx{i} = symfun(sym(['c_hx' suffix]), [x,y]) == c_hx_expr;
                
                c_0x_expr(x,y) = (1/obj.tau_x) * int((1/obj.tau_y) * obj.B_funcs{i}, x, obj.xl, obj.xl + obj.tau_x);
                obj.eq_e_ny = symfun(sym(['c_0x' suffix]), [x,y]) == c_0x_expr;
            end
            
            for i = 1:length(obj.T_funcs)
                d_hx_expr(x,y) = (2/obj.tau_x) * int(obj.T_funcs{i} .* cos(obj.beta_h*(x - obj.xl)), x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
                obj.eq_d_hx{i} = symfun(sym(['d_hx' suffix]), [x,y]) == d_hx_expr;
                
                d_0x_expr(x,y) = (1/obj.tau_x) * int((1/obj.tau_y) * obj.T_funcs{i}, x, obj.xl, obj.xl + obj.tau_x);
                obj.eq_f_ny{i} = symfun(sym(['d_0x' suffix]), [x,y]) == d_0x_expr;
            end
            
            for i = 1:length(obj.L_funcs)
                e_ny_expr(x,y) = (2/obj.tau_y) * int(obj.L_funcs{i} .* sin(obj.lambda_n*(y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y,'Hold',true);
                obj.eq_e_ny{i} = symfun(sym(['e_ny' suffix]), [x,y]) == e_ny_expr;
            end
            
            for i = 1:length(obj.R_funcs)
                f_ny_expr(x,y) = (2/obj.tau_y) * int(obj.R_funcs{i} .* sin(obj.lambda_n*(y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y,'Hold',true);
                obj.eq_f_ny{i} = symfun(sym(['f_ny' suffix]), [x,y]) == f_ny_expr;
            end
            
        end
    end
end
