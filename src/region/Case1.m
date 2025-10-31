classdef Case1 < BasicCase
    methods
        function obj = Case1(idx, xl, xr, yl, yt)
            % 调用父类构造函数
            obj@BasicCase(idx, xl, xr, yl, yt);
        end
        
        function regions = gen_solution_func(obj, Ln, Rn, Tn, Bn)
            syms x y
            % 应用边界条件
            obj.apply_boundaries(Ln, Rn, Tn, Bn);
            % 构造符号求和
            % 沿 x 方向的 A_z 分量
            
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
            
            % 总的 A_z
            A_z = A_zx_expr + A_zy_expr;
            
            % 计算磁场分量
            B_x = diff(A_z, y);
            B_y = -diff(A_z, x);
            
            % 构造符号方程
            suffix = ['_' num2str(obj.idx)];
            eq_A_z = symfun(sym(['A_z' suffix]), [x y]) == A_z;
            eq_B_x = symfun(sym(['B_x' suffix]), [x y]) == B_x + B_y;
            
            % 返回结果列表
            regions = {eq_A_z, eq_B_x, obj.eq_c_hx, obj.eq_d_hx, obj.eq_e_ny, obj.eq_f_ny};
        end
        
        % 设置系数方程
        function regions = gen_coefficient_func(obj)
            syms x y
            % 首先需要知道传入的边界方程，然后直接带入这个方程，这样子就是直接将边界方程带入下面的等式中
            % 也就是经过这个计算之后就直接确定了最后的参数方程
            % 计算 c_hx, d_hx, e_ny, f_ny 的积分表达式
            c_hx_expr(x,y) = (2 / obj.tau_x) * int(obj.B * sin(obj.beta_h * (x - obj.xl)), x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
            d_hx_expr(x,y) = (2 / obj.tau_x) * int(obj.T * sin(obj.beta_h * (x - obj.xl)), x, obj.xl, obj.xl + obj.tau_x,'Hold',true);
            e_ny_expr(x,y) = (2 / obj.tau_y) * int(obj.R * sin(obj.lambda_n * (y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y,'Hold',true);
            f_ny_expr(x,y) = (2 / obj.tau_y) * int(obj.L * sin(obj.lambda_n * (y - obj.yl)), y, obj.yl, obj.yl + obj.tau_y,'Hold',true);
            
            % 符号方程封装
            suffix = ['_' num2str(obj.idx)];
            obj.eq_c_hx = symfun(sym(['d_hx' suffix]), [x,y]) == c_hx_expr;
            obj.eq_d_hx = symfun(sym(['d_hx' suffix]), [x,y]) == d_hx_expr;
            obj.eq_e_ny = symfun(sym(['e_ny' suffix]), [x,y]) == e_ny_expr;
            obj.eq_f_ny = symfun(sym(['f_ny' suffix]), [x,y]) == f_ny_expr;
            
            regions = {obj.eq_c_hx, obj.eq_d_hx, obj.eq_e_ny, obj.eq_f_ny};
        end
    end
end
