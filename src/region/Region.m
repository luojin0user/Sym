classdef Region < handle
    
    properties
        idx
        mu_r
        J_r
        case_type       % CaseType 枚举类型
        boundary_status % 存储边界类型
        impl            % BasicCase 子类实例
        boundarys       % 边界管理对象
        all_regions     % 储存所有区域的节点，用于处理边界情况
        Ln
        Rn
        Tn
        Bn
    end
    
    methods
        function obj = Region(idx, type, xl, xr, yl, yt, bc_type, top, bottom, left, right, H_max, N_max, all_regions_num)
            % 修改可见边信息
            % L/R/T/B 为 bool 值，表示边界是否存在，如果为True，说明这个边界不存在
            obj.Ln = left(1)  == 0;
            obj.Rn = right(1) == 0;
            obj.Tn = top(1)   == 0;
            obj.Bn = bottom(1)== 0;
            
            obj.case_type = type;
            obj.idx = idx;
            
            
            switch type
                case CaseType.AlleyAir
                    obj.impl = AlleyAir(idx, xl, xr, yl, yt, obj.Ln, obj.Rn, obj.Tn, obj.Bn, H_max, N_max);
                case CaseType.BTAir % NormalAir
                    obj.impl = BTAir(idx, xl, xr, yl, yt, obj.Ln, obj.Rn, obj.Tn, obj.Bn, H_max, N_max);
                case CaseType.NormalAir
                    obj.impl = NormalAir(idx, xl, xr, yl, yt, obj.Ln, obj.Rn, obj.Tn, obj.Bn, H_max, N_max);
                case CaseType.FerriteCurrent
                    obj.impl = FerriteCurrent(idx, xl, xr, yl, yt, obj.Ln, obj.Rn, obj.Tn, obj.Bn, H_max, N_max);
                case CaseType.Aluminum
                    obj.impl = Aluminum(idx, xl, xr, yl, yt, obj.Ln, obj.Rn, obj.Tn, obj.Bn, H_max, N_max);
                otherwise
                    error('Unsupported case type');
            end
            
            obj.impl.num_idx_hn = zeros(all_regions_num, 2, 'uint32');
            
            obj.boundarys = Boundarys(obj, all_regions_num);
            obj.set_boundary(top, bottom, left, right, bc_type);
        end
        
        % 设置边界邻接的区域，每个传入一个数组
        function set_boundary(obj, top, bottom, left, right, bc_type)
            % 手动设置边界邻接信息
            obj.boundarys.bottom = bottom;
            obj.boundarys.top = top;
            obj.boundarys.left = left;
            obj.boundarys.right = right;
            obj.boundarys.bc_type = bc_type;
            
            % 同时将对应的边写入BasicCase类中
            obj.impl.tops = top;
            obj.impl.bottoms = bottom;
            obj.impl.lefts = left;
            obj.impl.rights = right;
            
        end
        
        
        % 计算这个区域内部的方程，这里直接返回这个方程，计算过程由set_boundary中的obj.impl.gen_region完成
        function region_func = get_region_solution_func(obj)
            obj.impl.gen_solution_func();
            region_func = {obj.impl.eq_A_z, obj.impl.eq_B_x};
        end
        
        
        % 计算这个区域系数的方程，这里是在所有区域内部方程计算完成之后调用的
        function funcs = gen_region_coefficient_func(obj, all_regions)
            obj.all_regions = all_regions;
            obj.impl.BCfuncs_loc_map = zeros(2, length(all_regions));
            obj.impl.all_regions = all_regions;
            % 首先需要计算边界方程，调用Cal_BC方法，然后就修改了Region中的impl属性的边界方程等等
            obj.boundarys.cal_BC(obj.Ln, obj.Rn, obj.Tn, obj.Bn);
            % 获得了边界方程之后，调用impl的方法生成
            obj.impl.gen_coefficient_func();
            % 收集系数方程及其对应的值，上下一一对应
            funcs = [obj.impl.eq_c0x; obj.impl.eq_c_hx; obj.impl.eq_d0x; obj.impl.eq_d_hx; obj.impl.eq_e_ny; obj.impl.eq_f_ny];
        end
        
    end
end
