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
                    eqT = obj.genAA(top_idx);
                    obj.impl.impl.T = eqT;
                    
                    % 下边界
                    bottom_idx = 0;
                    if Bn ~= 0 && ~isempty(obj.bottom)
                        bottom_idx = obj.bottom(1);
                    end
                    eqB = obj.genAA(bottom_idx);
                    obj.impl.impl.B = eqB;
                    
                    % 左边界
                    left_idx = 0;
                    if Ln ~= 0 && ~isempty(obj.left)
                        left_idx = obj.left(1);
                    end
                    eqL = obj.genAA(left_idx);
                    obj.impl.impl.L = eqL;
                    
                    % 右边界
                    right_idx = 0;
                    if Rn ~= 0 && ~isempty(obj.right)
                        right_idx = obj.right(1);
                    end
                    eqR = obj.genAA(right_idx);
                    obj.impl.impl.R = eqR;
                    
                case BC_TYPE.AABB
                    % 暂时空
            end
            
            bc_eqs = {eqT, eqB, eqL, eqR};
        end
        
        function eqF = genAA(obj, right_idx)
            global x y
            suffix = ['_' num2str(obj.impl.idx)];
            F = symfun(sym(['A_z' suffix]), [x y]);
            if right_idx == 0
                eqF = 0;
            else
                % 邻接区域的 A_z
                % 当是邻接区域的，直接把这个邻接区域的方程A_z送给对应的top
                G = obj.impl.all_regions{right_idx}.region_func{1};     % region_func的第1个就是对应的Az
                eqF = rhs(G);
            end
        end
    end
end
