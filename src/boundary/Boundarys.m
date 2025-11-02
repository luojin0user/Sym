classdef Boundarys < handle
    properties
        top
        bottom
        left
        right
        bc_type = BC_TYPE.AAAA
        impl
    end
    
    methods
        function obj = Boundarys(impl)
            obj.impl = impl;
            
        end
        
        function set_top(obj, l)
            obj.top = l;
        end
        
        function set_bottom(obj, l)
            obj.bottom = l;
        end
        
        function set_left(obj, l)
            obj.left = l;
        end
        
        function set_right(obj, l)
            obj.right = l;
        end
        
        function bc_eqs = cal_BC(obj, Ln, Rn, Tn, Bn)
            syms x y real
            % 根据 bc_type 生成边界方程
            switch obj.bc_type
                case BC_TYPE.BBAA
                    % 暂时空
                case BC_TYPE.AAAA
                    % 上边界
                    top_idx = 0;
                    if Tn ~= 0 && ~isempty(obj.top)
                        top_idx = obj.top(1);
                    end
                    [coeffs, eqT] = obj.genAA(top_idx);
                    for i = 1:length(eqT)
                        eqT{i} = subs(eqT{i}, y, obj.impl.impl.yt);
                    end
                    obj.impl.impl.T_funcs = eqT;
                    obj.impl.impl.T_coeffs = coeffs;
                    
                    % 下边界
                    bottom_idx = 0;
                    if Bn ~= 0 && ~isempty(obj.bottom)
                        bottom_idx = obj.bottom(1);
                    end
                    [coeffs, eqB] = obj.genAA(bottom_idx);
                    for i = 1:length(eqB)
                        eqB{i} = subs(eqB{i}, y, obj.impl.impl.yl);
                    end
                    obj.impl.impl.B_funcs = eqB;
                    obj.impl.impl.B_coeffs = coeffs;
                    
                    
                    % 左边界
                    left_idx = 0;
                    if Ln ~= 0 && ~isempty(obj.left)
                        left_idx = obj.left(1);
                    end
                    [coeffs, eqL] = obj.genAA(left_idx);
                    for i = 1:length(eqL)
                        eqL{i} = subs(eqL{i}, x, obj.impl.impl.xl);
                    end
                    obj.impl.impl.L_funcs = eqL;
                    obj.impl.impl.L_coeffs = coeffs;
                    
                    
                    % 右边界
                    right_idx = 0;
                    if Rn ~= 0 && ~isempty(obj.right)
                        right_idx = obj.right(1);
                    end
                    [coeffs, eqR] = obj.genAA(right_idx);
                    for i = 1:length(eqR)
                        eqR{i} = subs(eqR{i}, x, obj.impl.impl.xr);
                    end
                    obj.impl.impl.R_funcs = eqR;
                    obj.impl.impl.R_coeffs = coeffs;
                    
                case BC_TYPE.AABB
                    % 暂时空
            end
            
            bc_eqs = {eqT, eqB, eqL, eqR};
        end
        
        function [coeffs, eqF] = genAA(obj, right_idx)
            % 首先需要找到对应的方程，如果对应的方程中间包含多个c d e f，需要一一进行判断处理，然后送入数组中，和这个区域的边界情况的方程匹配
            % 如果由多个，送入数组中，由region类对这个数组进行处理
            eqF = {};
            coeffs = {};
            if right_idx == 0
                eqF = {0};
            else
                % 邻接区域的 A_z
                % 当是邻接区域的，直接把这个邻接区域的方程A_z送给对应的top
                % 除此之外，还要令对应的方程的变量值为边界值
                edge_region = obj.impl.all_regions{right_idx}; % 获取邻接区域对象
                edge_impl = edge_region.impl;
                if edge_region.B ~= 0
                    eqF{end+1} = subs(edge_impl.A_zx_expr, edge_impl.d_hx, 0);
                    eqF{end} = subs(edge_impl.A_zx_expr, edge_impl.c_hx, 1);
                    coeffs{end+1} = edge_impl.c_hx;
                end
                if edge_region.T ~= 0
                    eqF{end+1} = subs(edge_impl.A_zx_expr, edge_impl.c_hx, 0);
                    eqF{end} = subs(edge_impl.A_zx_expr, edge_impl.d_hx, 1);
                    coeffs{end+1} = edge_impl.d_hx;
                end
                if edge_region.R ~= 0
                    eqF{end+1} = subs(edge_impl.A_zy_expr, edge_impl.f_ny, 0);
                    eqF{end} = subs(edge_impl.A_zy_expr, edge_impl.e_ny, 1);
                    coeffs{end+1} = edge_impl.e_ny;
                end
                if edge_region.L ~= 0
                    eqF{end+1} = subs(edge_impl.A_zy_expr, edge_impl.e_ny, 0);
                    eqF{end} = subs(edge_impl.A_zy_expr, edge_impl.f_ny, 1);
                    coeffs{end+1} = edge_impl.f_ny;
                end
                
                % G = obj.impl.all_regions{right_idx}.region_func{1};     % region_func的第1个就是对应的Az
            end
        end
    end
end
