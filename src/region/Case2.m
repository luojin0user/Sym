classdef Case2 < BasicCase
    properties
        J_z
        c_0x
        d_0x
        A_y_P
    end
    
    methods
        function obj = Case2(idx, xl, xr, yl, yt)
            % 调用父类构造函数
            obj@BasicCase(idx, xl, xr, yl, yt);
            
            suffix = ['_' num2str(idx)];
            syms(['J_z' suffix], ['c_0x' suffix], ['d_0x' suffix], 'real')
            
            obj.J_z = eval(['J_z' suffix]);
            obj.c_0x = eval(['c_0x' suffix]);
            obj.d_0x = eval(['d_0x' suffix]);
            
            % 线性项
            obj.A_y_P = 0.5 * obj.mu_0 * obj.mu_r * obj.J_z * y^2;
        end
        
        function regions = gen_solution_func(obj, Ln, Rn, Tn, Bn)
            syms x y
            % 应用边界条件
            obj.apply_boundaries(Ln, Rn, Tn, Bn);
            
            % 构造符号求和
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
            
            A_z = A_zx_expr + A_zy_expr;
            
            % 磁场分量
            B_x = diff(A_z, y);
            B_y = -diff(A_z, x);
            
            % 符号方程
            suffix = ['_' num2str(obj.idx)];
            eq_A_z = A_z == symfun(sym(['A_z' suffix]), [x y]);
            eq_B_x = (B_x + B_y) == symfun(sym(['B_x' suffix]), [x y]);
            
            % 返回 cell array
            regions = {eq_A_z, eq_B_x, c_hx_int, d_hx_int, e_ny_int, f_ny_int, c_0x_int, d_0x_int};
        end
        
        
        % 设置系数方程
        function regions = gen_coefficient_func(obj)
            syms x y
            % 傅里叶系数积分
            c_0x_expr(x,y) = (1/obj.tau_x) * int((1/obj.tau_y) * obj.B, x, obj.xl, obj.xl + obj.tau_x);
            d_0x_expr(x,y) = (1/obj.tau_x) * int((1/obj.tau_y) * obj.T, x, obj.xl, obj.xl + obj.tau_x);
            c_hx_expr(x,y) = (2/obj.tau_x) * int(obj.B .* cos(obj.beta_h*(x - obj.xl)), x, obj.xl, obj.xl + obj.tau_x);
            d_hx_expr(x,y) = (2/obj.tau_x) * int(obj.T .* cos(obj.beta_h*(x - obj.xl)), x, obj.xl, obj.xl + obj.tau_x);
            e_ny_expr(x,y) = (2/obj.tau_y) * int(obj.R .* sin(obj.lambda_n*(y - yl)), y, yl, yl + obj.tau_y);
            f_ny_expr(x,y) = (2/obj.tau_y) * int(obj.L .* sin(obj.lambda_n*(y - yl)), y, yl, yl + obj.tau_y);
            
            % 符号方程封装
            suffix = ['_' num2str(obj.idx)];
            obj.eq_e_ny = symfun(sym(['c_0x' suffix]), [x,y]) == c_0x_expr;
            obj.eq_f_ny = symfun(sym(['d_0x' suffix]), [x,y]) == d_0x_expr;
            
            obj.eq_c_hx = symfun(sym(['d_hx' suffix]), [x,y]) == c_hx_expr;
            obj.eq_d_hx = symfun(sym(['d_hx' suffix]), [x,y]) == d_hx_expr;
            obj.eq_e_ny = symfun(sym(['e_ny' suffix]), [x,y]) == e_ny_expr;
            obj.eq_f_ny = symfun(sym(['f_ny' suffix]), [x,y]) == f_ny_expr;
            
            regions = {obj.eq_c_hx, obj.eq_d_hx, obj.eq_e_ny, obj.eq_f_ny};
        end
    end
end
