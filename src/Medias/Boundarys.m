classdef Boundarys < handle
    properties
        top
        bottom
        left
        right
        bc_type
        impl
        case_impl
    end
    
    methods
        function obj = Boundarys(impl)
            obj.impl = impl;
            obj.case_impl = impl.impl;  % 这个指的是如BTAir类的一个实例
        end
        
        
        function cal_BC(obj, Ln, Rn, Tn, Bn)
            syms x y real
            
            % 根据 bc_type 生成边界方程
            
            switch obj.bc_type
                case BC_TYPE.BBAA
                    % 这里需要注意可能出现的分段函数的形式，如果出现分段函数，一定是B连续，所以需要
                    % 目前而言，出现分段函数只会在上下边界为BB的情况出现
                    % 当出现分段时，采用的储存方法为BC{i}{j}，BC{1}{j}是上边的第一个分段区域，BC{2}{j}是上边的第二个分段区域，以此类推
                    % 这样表示的方法和论文中的类似
                    % 上边界
                    if Tn == 0 && ~isempty(obj.top)
                        for i = 1:length(obj.top)
                            % 每次取一个，计算他的边界函数
                            top_idx = obj.top(i);
                            [coeffs, eqT] = obj.genBB(top_idx, 1); % 这里取得的eqT是一个cell数组，其中的每一项是边界方程的一部分，例如c这部分，
                            for j = 1:length(eqT)
                                eqT{j} = subs(eqT{j}, y, obj.case_impl.yt);
                            end
                            % 如果是分段函数，传入的T_funcs将会是一个更长的一维cell数组
                            obj.case_impl.T_funcs = [obj.case_impl.T_funcs, eqT];     % 这样子就成为了一个一维数组，每次传入其中的一行
                            obj.case_impl.T_coeffs = [obj.case_impl.T_coeffs ,coeffs];
                        end
                    else
                        obj.case_impl.T_funcs = {};
                        obj.case_impl.T_coeffs = {};
                    end
                    
                    % 下边界
                    if Bn == 0 && ~isempty(obj.bottom)
                        for i = 1:length(obj.bottom)
                            % 每次取一个，计算他的边界函数
                            bottom_idx = obj.bottom(i);
                            [coeffs, eqB] = obj.genBB(bottom_idx, 1); % 这里取得的eqT是一个cell数组，其中的每一项是边界方程的一部分，例如c这部分，
                            for j = 1:length(eqB)
                                eqB{j} = subs(eqB{j}, y, obj.case_impl.yl);
                            end
                            % 如果是分段函数，传入的T_funcs将会是一个更长的一维cell数组
                            obj.case_impl.B_funcs = [obj.case_impl.B_funcs, eqB];     % 这样子就成为了一个一维数组，每次传入其中的一行
                            obj.case_impl.B_coeffs = [obj.case_impl.B_coeffs ,coeffs];
                        end
                    else
                        obj.case_impl.B_funcs = {};
                        obj.case_impl.B_coeffs = {};
                    end
                    
                    % 左边界
                    if Ln == 0 && ~isempty(obj.left)
                        left_idx = obj.left(1);
                        
                        [coeffs, eqL] = obj.genAA(left_idx);
                        for i = 1:length(eqL)
                            eqL{i} = eqL{i} * obj.case_impl.lambda_n;
                            eqL{i} = subs(eqL{i}, x, obj.case_impl.xl);
                        end
                        obj.case_impl.L_funcs = eqL;
                        obj.case_impl.L_coeffs = coeffs;
                    else
                        obj.case_impl.L_funcs = {};
                        obj.case_impl.L_coeffs = {};
                    end
                    
                    % 右边界
                    if Rn == 0 && ~isempty(obj.right)
                        right_idx = obj.right(1);
                        
                        [coeffs, eqR] = obj.genAA(right_idx);
                        for i = 1:length(eqR)
                            eqR{i} = eqR{i} * obj.case_impl.lambda_n;
                            eqR{i} = subs(eqR{i}, x, obj.case_impl.xr);
                        end
                        obj.case_impl.R_funcs = eqR;
                        obj.case_impl.R_coeffs = coeffs;
                    else
                        obj.case_impl.R_funcs = {};
                        obj.case_impl.R_coeffs = {};
                    end
                    
                case BC_TYPE.AAAA
                    % 上边界
                    if Tn == 0 && ~isempty(obj.top)
                        top_idx = obj.top(1);
                        [coeffs, eqT] = obj.genAA(top_idx);
                        for i = 1:length(eqT)
                            eqT{i} = eqT{i} * obj.case_impl.beta_h;  % 乘以 beta_h
                            eqT{i} = subs(eqT{i}, y, obj.case_impl.yt);
                        end
                        obj.case_impl.T_funcs = eqT;
                        obj.case_impl.T_coeffs = coeffs;
                    else
                        obj.case_impl.T_funcs = {};
                        obj.case_impl.T_coeffs = {};
                    end
                    
                    % 下边界
                    if Bn == 0 && ~isempty(obj.bottom)
                        bottom_idx = obj.bottom(1);
                        
                        [coeffs, eqB] = obj.genAA(bottom_idx);
                        for i = 1:length(eqB)
                            eqB{i} = eqB{i} * obj.case_impl.beta_h;  % 乘以 beta_h
                            eqB{i} = subs(eqB{i}, y, obj.case_impl.yl);
                        end
                        obj.case_impl.B_funcs = eqB;
                        obj.case_impl.B_coeffs = coeffs;
                    else
                        obj.case_impl.B_funcs = {};
                        obj.case_impl.B_coeffs = {};
                    end
                    
                    % 左边界
                    if Ln == 0 && ~isempty(obj.left)
                        left_idx = obj.left(1);
                        
                        [coeffs, eqL] = obj.genAA(left_idx);
                        for i = 1:length(eqL)
                            eqL{i} = eqL{i} * obj.case_impl.lambda_n;
                            eqL{i} = subs(eqL{i}, x, obj.case_impl.xl);
                        end
                        obj.case_impl.L_funcs = eqL;
                        obj.case_impl.L_coeffs = coeffs;
                    else
                        obj.case_impl.L_funcs = {};
                        obj.case_impl.L_coeffs = {};
                    end
                    
                    % 右边界
                    if Rn == 0 && ~isempty(obj.right)
                        right_idx = obj.right(1);
                        
                        [coeffs, eqR] = obj.genAA(right_idx);
                        for i = 1:length(eqR)
                            eqR{i} = eqR{i} * obj.case_impl.lambda_n;
                            eqR{i} = subs(eqR{i}, x, obj.case_impl.xr);
                        end
                        obj.case_impl.R_funcs = eqR;
                        obj.case_impl.R_coeffs = coeffs;
                    else
                        obj.case_impl.R_funcs = {};
                        obj.case_impl.R_coeffs = {};
                    end
                    
                case BC_TYPE.AABB
                    % 这里只考虑了非分段函数的情况，即传入(1)，后续需要考虑分段函数，传入所有的邻接边
                    % 上边界
                    if Tn == 0 && ~isempty(obj.top)
                        top_idx = obj.top(1);
                        [coeffs, eqT] = obj.genAA(top_idx);
                        for i = 1:length(eqT)
                            eqT{i} = subs(eqT{i}, y, obj.case_impl.yt);
                        end
                        obj.case_impl.T_funcs = eqT;
                        obj.case_impl.T_coeffs = coeffs;
                    else
                        obj.case_impl.T_funcs = {};
                        obj.case_impl.T_coeffs = {};
                    end
                    
                    % 下边界
                    if Bn == 0 && ~isempty(obj.bottom)
                        bottom_idx = obj.bottom(1);
                        
                        [coeffs, eqB] = obj.genAA(bottom_idx);
                        for i = 1:length(eqB)
                            eqB{i} = subs(eqB{i}, y, obj.case_impl.yl);
                        end
                        obj.case_impl.B_funcs = eqB;
                        obj.case_impl.B_coeffs = coeffs;
                    else
                        obj.case_impl.B_funcs = {};
                        obj.case_impl.B_coeffs = {};
                    end
                    
                    % 左边界
                    if Ln == 0 && ~isempty(obj.left)
                        left_idx = obj.left(1);
                        
                        [coeffs, eqL] = obj.genBB(left_idx, 2); % 左右边界是ef系数
                        for i = 1:length(eqL)
                            eqL{i} = subs(eqL{i}, x, obj.case_impl.xl);
                        end
                        obj.case_impl.L_funcs = eqL;
                        obj.case_impl.L_coeffs = coeffs;
                    else
                        obj.case_impl.L_funcs = {};
                        obj.case_impl.L_coeffs = {};
                    end
                    
                    % 右边界
                    if Rn == 0 && ~isempty(obj.right)
                        right_idx = obj.right(1);
                        
                        [coeffs, eqR] = obj.genBB(right_idx, 2);
                        for i = 1:length(eqR)
                            eqR{i} = subs(eqR{i}, x, obj.case_impl.xr);
                        end
                        obj.case_impl.R_funcs = eqR;
                        obj.case_impl.R_coeffs = coeffs;
                    else
                        obj.case_impl.R_funcs = {};
                        obj.case_impl.R_coeffs = {};
                    end
            end
            % 这里修改了impl中的属性，无需再返回值
        end
        
        function [coeffs, eqF] = genAA(obj, right_idx)
            % 首先需要找到对应的方程，如果对应的方程中间包含多个c d e f，需要一一进行判断处理，然后送入数组中，和这个区域的边界情况的方程匹配
            % 如果由多个，送入数组中，由region类对这个数组进行处理
            eqF = cell(1,6);
            coeffs = cell(1,6);
            if right_idx == 0
                eqF = {0};
            else
                % 邻接区域的 A_z，
                % 当是邻接区域的，直接把这个邻接区域的方程A_z送给对应的top
                % 除此之外，还要令对应的方程的变量值为边界值
                edge_region = obj.impl.all_regions{right_idx}; % 获取邻接区域对象
                edge_impl = edge_region.impl;
                has_cd0x = (edge_region.case_type == CaseType.FerriteCurrent);     % 当邻接区域是FerriteCurrent时，需要处理
                if has_cd0x
                    % 首先处理c_0x
                    F = subs(edge_impl.A_zx_expr, {edge_impl.c_hx, edge_impl.d_hx, edge_impl.c_0x, edge_impl.d_0x}, {0, 0, 1, 0});    % 只留下c_0x
                    eqF{5} = F;
                    coeffs{5} = edge_impl.c_0x;
                    
                    F = subs(edge_impl.A_zx_expr, {edge_impl.c_hx, edge_impl.d_hx, edge_impl.c_0x, edge_impl.d_0x}, {0, 0, 0, 1});    % 只留下c_0x
                    eqF{6} = F;
                    coeffs{6} = edge_impl.d_0x;
                end
                
                if edge_region.Bn == 0  % 相邻区域的这个边界上有c_hx
                    % 如果有c_0x，说明这个需要考虑分段函数
                    % 分段函数直接放在最后两个位置，分别是c_0x和d_0x的表达式
                    F1 = subs(edge_impl.A_zx_expr, {edge_impl.c_hx, edge_impl.d_hx}, {1, 0});
                    if has_cd0x
                        F1 = subs(F1, {edge_impl.c_0x, edge_impl.d_0x}, {0, 0});    % 如果有c_0x,d_0x的话，置为0
                    end
                    eqF{1} = F1;
                    coeffs{1} = edge_impl.c_hx;
                end
                if edge_region.Tn == 0
                    F1 = subs(edge_impl.A_zx_expr, {edge_impl.c_hx, edge_impl.d_hx}, {0, 1});
                    if has_cd0x
                        F1 = subs(F1, {edge_impl.c_0x, edge_impl.d_0x}, {0, 0});    % 如果有c_0x,d_0x的话，置为0
                    end
                    eqF{2} = F1;
                    coeffs{2} = edge_impl.d_hx;
                end
                if edge_region.Ln == 0
                    F1 = subs(edge_impl.A_zy_expr, {edge_impl.e_ny, edge_impl.f_ny}, {1, 0});
                    eqF{3} = F1;
                    coeffs{3} = edge_impl.e_ny;
                end
                if edge_region.Rn == 0
                    F1 = subs(edge_impl.A_zy_expr, {edge_impl.e_ny, edge_impl.f_ny}, {0, 1});
                    eqF{4} = F1;
                    coeffs{4} = edge_impl.f_ny;
                end
                
                
                % G = obj.impl.all_regions{right_idx}.region_func{1};     % region_func的第1个就是对应的Az
            end
        end
        
        
        function [coeffs, eqF] = genBB(obj, right_idx, cd_or_ef)
            % 这里记得乘以2个的mu_0的系数
            % 首先需要找到对应的方程，如果对应的方程中间包含多个c d e f，需要一一进行判断处理，然后送入数组中，和这个区域的边界情况的方程匹配
            % 如果由多个，送入数组中，由region类对这个数组进行处理
            eqF = cell(1,6);
            coeffs = cell(1,6);
            if right_idx == 0
                eqF = {0};
            else
                % 邻接区域的 A_z，
                % 当是邻接区域的，直接把这个邻接区域的方程A_z送给对应的top
                % 除此之外，还要令对应的方程的变量值为边界值
                % 对于一个已经确定的区域，其e、f参数由B_y决定，其中某一个(e)的参数c d由B_y_x决定，e f由B_y_y决定
                edge_region = obj.impl.all_regions{right_idx}; % 获取邻接区域对象
                edge_impl = edge_region.impl;
                has_cd0x = (edge_region.case_type == CaseType.FerriteCurrent);     % 当邻接区域是FerriteCurrent时，需要处理
                if has_cd0x
                    F1 = obj.getBB_func(edge_impl, 1, cd_or_ef);
                    % 首先处理c_0x
                    F = subs(F1, {edge_impl.c_hx, edge_impl.d_hx, edge_impl.c_0x, edge_impl.d_0x}, {0, 0, 1, 0});    % 只留下c_0x
                    eqF{5} = F;
                    coeffs{5} = edge_impl.c_0x;
                    
                    F = subs(F1, {edge_impl.c_hx, edge_impl.d_hx, edge_impl.c_0x, edge_impl.d_0x}, {0, 0, 0, 1});    % 只留下c_0x
                    eqF{6} = F;
                    coeffs{6} = edge_impl.d_0x;
                end
                if edge_region.Bn == 0  % 相邻区域的这个边界上有c_hx
                    % 如果有c_0x，说明这个需要考虑分段函数
                    % 分段函数直接放在最后两个位置，分别是c_0x和d_0x的表达式
                    F1 = obj.getBB_func(edge_impl, 1, cd_or_ef);
                    F1 = subs(F1, {edge_impl.c_hx, edge_impl.d_hx}, {1, 0});
                    if has_cd0x
                        F1 = subs(F1, {edge_impl.c_0x, edge_impl.d_0x}, {0, 0});    % 如果有c_0x,d_0x的话，置为0
                    end
                    eqF{1} = F1;
                    coeffs{1} = edge_impl.c_hx;
                end
                if edge_region.Tn == 0  % 相邻区域的这个边界上有d_hx
                    % 如果有c_0x，说明这个需要考虑分段函数
                    % 分段函数直接放在最后两个位置，分别是c_0x和d_0x的表达式
                    F1 = obj.getBB_func(edge_impl, 1, cd_or_ef);
                    F1 = subs(F1, {edge_impl.c_hx, edge_impl.d_hx}, {0, 1});
                    if has_cd0x
                        F1 = subs(F1, {edge_impl.c_0x, edge_impl.d_0x}, {0, 0});    % 如果有c_0x,d_0x的话，置为0
                    end
                    eqF{2} = F1;
                    coeffs{2} = edge_impl.d_hx;
                end
                if edge_region.Ln == 0
                    F1 = obj.getBB_func(edge_impl, 2, cd_or_ef);
                    F = subs(F1, {edge_impl.e_ny, edge_impl.f_ny}, {1,0});
                    eqF{3} = F;
                    coeffs{3} = edge_impl.e_ny;
                end
                
                if edge_region.Rn == 0
                    F1 = obj.getBB_func(edge_impl, 2, cd_or_ef);
                    F = subs(F1, {edge_impl.e_ny, edge_impl.f_ny}, {0,1});
                    eqF{4} = F;
                    coeffs{4} = edge_impl.f_ny;
                end
                
                
                for i = 1:length(eqF)
                    if ~isempty(eqF{i})
                        eqF{i} = eqF{i} * (obj.case_impl.mu_r / edge_impl.mu_r);  % 乘以 mu_r / mu_r
                    end
                end
                
                % G = obj.impl.all_regions{right_idx}.region_func{1};     % region_func的第1个就是对应的Az
            end
        end
        
        function F1 = getBB_func(obj, edge_impl, x_or_y, cd_or_ef)
            if( x_or_y == 1)
                if cd_or_ef == 1
                    F1 = edge_impl.B_x_x;
                else
                    F1 = edge_impl.B_y_x;
                end
            else
                if cd_or_ef == 1
                    F1 = edge_impl.B_x_y;
                else
                    F1 = edge_impl.B_y_y;
                end
            end
        end
    end
end
