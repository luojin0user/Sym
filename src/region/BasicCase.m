classdef BasicCase < handle
    properties
        idx
        
        % --- 公共编号变量 ---
        xr
        xl
        yt
        yl
        tau_x
        tau_y
        
        
        % --- 公共函数 ---
        T
        B
        L
        R
        
        % --- 公共求和指标 ---
        h
        n
        
        % --- 公共参数 ---
        beta_h
        lambda_n
        mu_r
        H_max
        N_max
        
        % --- 边界积分项 ---
        c_hx
        d_hx
        e_ny
        f_ny
        
        eq_c0x
        eq_d0x
        eq_c_hx
        eq_d_hx
        eq_e_ny
        eq_f_ny
    end
    
    methods
        % 构造函数
        function obj = BasicCase(idx, xl, xr, yl, yt)
            global x y mu_0
            if nargin < 1
                error('必须指定 idx');
            end
            obj.idx = idx;
            suffix = ['_' num2str(idx)];
            
            % --- 公共编号变量 ---
            % 定义符号变量（不带后缀）
            %{
            xr = sym('xr', 'real');
            xl = sym('xl', 'real');
            yt = sym('yt', 'real');
            yl = sym('yl', 'real');
            %}
            obj.xr = xr;
            obj.xl = xl;
            obj.yt = yt;
            obj.yl = yl;
            
            obj.tau_x = obj.xr - obj.xl;
            obj.tau_y = obj.yt - obj.yl;
            
            % --- 公共函数 ---
            obj.T = symfun(sym('T'), [x, y]);  % T(x,y)
            obj.B = symfun(sym('B'), [x, y]);  % B(x,y)
            obj.L = symfun(sym('L'), [x, y]);  % L(x,y)
            obj.R = symfun(sym('R'), [x, y]);  % R(x,y)
            
            
            % --- 公共求和指标 ---
            syms h n integer;
            obj.h = h;
            obj.n = n;
            
            % 定义符号变量
            syms mu_r 'real'
            syms H_max N_max 'integer'
            
            % 定义符号函数
            obj.beta_h = h * pi / obj.tau_x;      % beta 是关于 h 的函数
            obj.lambda_n = n * pi / obj.tau_y;  % lambda 是关于 n 的函数
            
            % 直接赋值给对象属性
            obj.mu_r = mu_r;
            obj.H_max = H_max;
            obj.N_max = N_max;
            
            
            % --- 边界积分项 ---
            obj.c_hx = symfun(sym(['c_hx' suffix]), x);  % c_hx_1(x)
            obj.d_hx = symfun(sym(['d_hx' suffix]), x);  % d_hx_1(x)
            obj.e_ny = symfun(sym(['e_ny' suffix]), y);  % e_ny_1(y)
            obj.f_ny = symfun(sym(['f_ny' suffix]), y);  % f_ny_1(y)
        end
        
        % 公共逻辑：处理 Ln/Rn/Tn/Bn
        function apply_boundaries(obj, Ln, Rn, Tn, Bn)
            if ~Ln
                obj.f_ny = 0;
            end
            if ~Rn
                obj.e_ny = 0;
            end
            if ~Tn
                obj.d_hx = 0;
            end
            if ~Bn
                obj.c_hx = 0;
            end
        end
        
        % 用于生成这个区域内部的函数
        function regions = gen_solution_func(obj, Ln, Rn, Tn, Bn)
            error('子类必须实现 gen_region 方法');
        end
        
        % 用于生成这个区域的边界的函数
        function regions = gen_coefficient_func(obj)
            error('子类必须实现 gen_region 方法');
        end
    end
end
