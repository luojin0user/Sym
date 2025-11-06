classdef BasicCase < handle
    properties
        x
        y
        
        idx
        
        % --- 公共编号变量 ---
        xr
        xl
        yt
        yl
        tau_x
        tau_y
        
        
        % --- 公共函数 ---
        T_funcs
        B_funcs
        L_funcs
        R_funcs
        
        T_coeffs
        B_coeffs
        L_coeffs
        R_coeffs
        
        Ln
        Rn
        Tn
        Bn
        
        num_idx_hn  % num_idx_hn(idx,1)用于储存与idx区域边界相接的变量的h个数（对方的c d），ef个数为2
        
        % --- 公共求和指标 ---
        h
        n
        
        % --- 公共参数 ---
        beta_h
        lambda_n
        mu_0
        mu_r
        H_max
        N_max
        
        % --- 边界积分项 ---
        c_0x
        d_0x
        
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
        
        A_zx_expr
        A_zy_expr
        
        eq_A_z
        eq_B_x
        
        B_x_x
        B_x_y
        B_y_x
        B_y_y
        
        % 边界邻接区域
        tops
        bottoms
        lefts
        rights
        
        
        all_regions     % 储存所有区域
        num_coeffs      % 这个区域的参数个数
        coeffs_exists      % 当前区域的参数情况，1x6的数组，代表c0 c d0 d e f
        BCfuncs_loc_map   % 边界方程的位置，索引代表区域，值1代表边界方程的类型（2=c=B，4=d=T，5=e=L，6=f=R），值2代表方程索引
        edge_idx_
    end
    
    methods
        % 构造函数
        function obj = BasicCase(idx, xl, xr, yl, yt, H_max, N_max)
            syms x y
            if nargin < 1
                error('必须指定 idx');
            end
            obj.idx = idx;
            obj.x = x;
            obj.y = y;
            
            suffix = ['_' num2str(idx)];
            
            % --- 公共编号变量 ---
            % 定义符号变量（不带后缀）
            
            obj.xr = sym(['xr' suffix], 'real');
            obj.xl = sym(['xl' suffix], 'real');
            obj.yt = sym(['yt' suffix], 'real');
            obj.yl = sym(['yl' suffix], 'real');
            %{
            %}
            obj.xr = xr;
            obj.xl = xl;
            obj.yt = yt;
            obj.yl = yl;
            
            obj.tau_x = obj.xr - obj.xl;
            obj.tau_y = obj.yt - obj.yl;
            
            % --- 公共求和指标 ---
            obj.h = sym(['h' suffix], 'integer');
            obj.n = sym(['n' suffix], 'integer');
            
            % obj.h = H_max;
            % obj.n = N_max;
            
            % 定义符号变量
            syms mu_r 'real'
            % syms H_max N_max 'integer'
            
            % 定义符号函数
            obj.beta_h = obj.h * pi / obj.tau_x;      % beta 是关于 h 的函数
            obj.lambda_n = obj.n * pi / obj.tau_y;  % lambda 是关于 n 的函数
            
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
            obj.Ln = Ln;
            obj.Rn = Rn;
            obj.Tn = Tn;
            obj.Bn = Bn;
            
            
            if Ln   % 没有左边界
                obj.e_ny = 0;
                obj.num_coeffs = obj.num_coeffs - 1;
            end
            if Rn
                obj.f_ny = 0;
                obj.num_coeffs = obj.num_coeffs - 1;
            end
            if Tn   % 没有上边界
                obj.d_hx = 0;
                obj.num_coeffs = obj.num_coeffs - 1;
            end
            if Bn
                obj.c_hx = 0;
                obj.num_coeffs = obj.num_coeffs - 1;
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
