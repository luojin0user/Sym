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
        function obj = Region(c1_or_c2, idx, xl, xr, yl, yt, all_regions)
            if c1_or_c2 == 1
                obj.impl = Case1(idx, xl, xr, yl, yt);
            else
                obj.impl = Case2(idx, xl, xr, yl, yt);
            end
            obj.idx = idx;
            obj.all_regions = all_regions;
            obj.boundarys = Boundarys(obj); % 假设已有 MATLAB 版 Boundarys 类
            
        end
        
        function set_region(obj, mu_r, J_r)
            obj.mu_r = mu_r;
            obj.J_r = J_r;
            
            % 更新 BasicCase 内的几何参数
            obj.impl.xl = x1;
            obj.impl.xr = x2;
            obj.impl.yl = y1;
            obj.impl.yt = y2;
            obj.impl.tau_x = x2 - x1;
            obj.impl.tau_y = y2 - y1;
        end
        
        function set_H_N_max(obj, Hmax, Nmax)
            % obj.impl.H_max = Hmax;
            % obj.impl.N_max = Nmax;
        end
        
        function set_boundary(obj, L, R, T, B)
            % L/R/T/B 为 bool 值，表示边界是否存在
            obj.L = L;
            obj.R = R;
            obj.T = T;
            obj.B = B;
            
            % 生成区域方程
            obj.region_func = obj.impl.gen_region(L, R, T, B);
            
            % 示例：手动设置边界邻接信息
            obj.boundarys.bottom = [1];
            obj.boundarys.top = [1];
            obj.boundarys.left = [0];
            obj.boundarys.right = [1];
            
        end
        
        % 计算这个区域内部的方程，这里直接返回这个方程，计算过程由set_boundary中的obj.impl.gen_region完成
        function region_func = get_region_func(obj)
            region_func = obj.region_func;
            % region_bc_func = obj.region_bc_func;
        end
        
        % 计算这个区域系数的方程，这里是在所有区域内部方程计算完成之后调用的
        function regions = gen_region_coefficient_func(obj)
            % 首先需要计算边界方程，调用Cal_BC方法
            obj.region_bc_func = obj.boundarys.cal_BC(obj.L, obj.R, obj.T, obj.B);
            % 获得了边界方程之后，调用impl的方法生成
            regions = obj.impl.gen_coefficient_func();
        end
        
        function cal_boundary(obj, L, R, T, B)
            global x y
            suffix = ['_' num2str(obj.idx)];
            
            if L
                obj.impl.L = symfun(sym(['A_z' suffix]), [x y]);
            end
            if R
                obj.impl.R = symfun(sym(['A_z' suffix]), [x y]);
            end
            if T
                obj.impl.T = symfun(sym(['A_z' suffix]), [x y]);
            end
            if B
                obj.impl.B = symfun(sym(['A_z' suffix]), [x y]);
            end
        end
    end
end
