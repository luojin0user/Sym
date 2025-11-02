classdef Region < handle
    
    properties
        idx
        mu_r
        J_r
        c1_or_c2        % 1 = Case1, 2 = Case2
        region_func     % 区域内的方程列表
        region_bc_func  % 区域边界的方程列表
        boundary_status % 存储边界类型
        impl            % BasicCase 子类实例
        boundarys       % 边界管理对象
        all_regions     % 储存所有区域的节点，用于处理边界情况
        L
        R
        T
        B
    end
    
    methods
        function obj = Region(c1_or_c2, idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max, all_regions)
            % L/R/T/B 为 bool 值，表示边界是否存在
            obj.L = Ln;
            obj.R = Rn;
            obj.T = Tn;
            obj.B = Bn;
            
            if c1_or_c2 == 1
                obj.impl = BTAir(idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max);
            else
                obj.impl = NormalAir(idx, xl, xr, yl, yt, Ln, Rn, Tn, Bn, H_max, N_max);
            end
            obj.idx = idx;
            obj.all_regions = all_regions;
            obj.boundarys = Boundarys(obj); % 假设已有 MATLAB 版 Boundarys 类
            
        end
        
        
        function set_H_N_max(obj, Hmax, Nmax)
            obj.impl.H_max = Hmax;
            obj.impl.N_max = Nmax;
        end
        
        
        function set_boundary(obj)
            
            % 生成区域方程
            obj.region_func = obj.impl.gen_solution_func();
            
            % 示例：手动设置边界邻接信息
            obj.boundarys.bottom = [1];
            obj.boundarys.top = [1];
            obj.boundarys.left = [0];
            obj.boundarys.right = [1];
            
        end
        
        
        % 计算这个区域内部的方程，这里直接返回这个方程，计算过程由set_boundary中的obj.impl.gen_region完成
        function region_func = get_region_solution_func(obj)
            region_func = obj.region_func;
            % region_bc_func = obj.region_bc_func;
        end
        
        
        % 计算这个区域系数的方程，这里是在所有区域内部方程计算完成之后调用的
        function [funcs, means] = gen_region_coefficient_func(obj)
            % 首先需要计算边界方程，调用Cal_BC方法，然后就修改了Region中的impl属性的边界方程等等
            obj.boundarys.cal_BC(obj.L, obj.R, obj.T, obj.B);
            % 获得了边界方程之后，调用impl的方法生成
            obj.impl.gen_coefficient_func();
            % 收集系数方程及其对应的值，上下一一对应
            funcs = {obj.impl.eq_c_hx, obj.impl.eq_d_hx, obj.impl.eq_e_ny, obj.impl.eq_f_ny};
            means = {obj.impl.B_coeffs, obj.impl.T_coeffs, obj.impl.R_coeffs, obj.impl.L_coeffs};
        end
        
    end
end
