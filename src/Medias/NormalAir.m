classdef NormalAir < Case1
    methods
        function obj = NormalAir(idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max)
            % 调用父类构造函数
            obj@Case1(idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max);
            syms x y real
            
            obj.Ax_c = (1 / obj.beta_h) * sinh(obj.beta_h * (obj.yt - y)) / sinh(obj.beta_h * obj.tau_y);
            obj.Ax_d = (1 / obj.beta_h) * sinh(obj.beta_h * (y - obj.yl)) / sinh(obj.beta_h * obj.tau_y);
            obj.Ax_e = (1 / obj.lambda_n) * sinh(obj.lambda_n * (obj.xr - x)) / sinh(obj.lambda_n * obj.tau_x);
            obj.Ax_f = (1 / obj.lambda_n) * sinh(obj.lambda_n * (x - obj.xl)) / sinh(obj.lambda_n * obj.tau_x);
            
            obj.A_zx_expr = symsum( ...
                ( (obj.c_hx * obj.Ax_c) + (obj.d_hx * obj.Ax_d ) )* sin(obj.beta_h * (x - obj.xl)), ...
                obj.h, 1, obj.H_max);
            
            obj.A_zy_expr = symsum( ...
                ( (obj.e_ny * obj.Ax_e) + (obj.f_ny * obj.Ax_f)) * sin(obj.lambda_n * (y - obj.yl)), ...
                obj.n, 1, obj.N_max);
            
        end
    end
end
