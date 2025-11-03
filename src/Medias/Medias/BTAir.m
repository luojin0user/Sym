classdef BTAir < Case1
    methods
        function obj = BTAir(idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max)
            % 调用父类构造函数
            obj@Case1(idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max);
            syms x y real
            %{
            
            obj.A_zx_expr = symfun( symsum( ...
                ( (-obj.c_hx / obj.beta_h) * sinh(obj.beta_h * (obj.yt - y)) / cosh(obj.beta_h * obj.tau_y) ...
                + (obj.d_hx / obj.beta_h) * sinh(obj.beta_h * (y - obj.yl)) / cosh(obj.beta_h * obj.tau_y) ) ...
                * sin(obj.beta_h * (x - obj.xl)), ...
                obj.h, 1, obj.H_max), [x,y]);
            
            obj.A_zy_expr = 0;
            %}
            
            obj.A_zx_expr = symfun(  ...
                ( (-obj.c_hx / obj.beta_h) * sinh(obj.beta_h * (obj.yt - y)) / cosh(obj.beta_h * obj.tau_y) ...
                + (obj.d_hx / obj.beta_h) * sinh(obj.beta_h * (y - obj.yl)) / cosh(obj.beta_h * obj.tau_y) ) ...
                * sin(obj.beta_h * (x - obj.xl)), ...
                [x,y]);
            
            obj.A_zy_expr = 0;
        end
    end
end
