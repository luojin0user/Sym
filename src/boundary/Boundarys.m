classdef Boundarys < handle
    properties
        top
        bottom
        left
        right
        bc_type
        impl
        case_impl
        
        % 这个区域邻接方程的位置
        BCfuncs_loc % 第一列代表当前这个邻接区域对应是上（1）下（2）还是左（3）右（4），第二列是方程位置
        BCfuncs_loc_map % 记录<idx,is_cal>记录邻接区域的idx与
        all_region_num % 所有的区域数量
        
        region2edge % 储存区域到边界的映射，索引代表区域，值代表边界，其中1 2 3 4 5 6分别代表c0 c d0 d e f
    end
    
    methods
        function obj = Boundarys(impl, all_region_num)
            obj.impl = impl;
            obj.case_impl = impl.impl;  % 这个指的是如BTAir类的一个实例
            
            obj.all_region_num = all_region_num;
            obj.BCfuncs_loc = zeros(2, all_region_num);     % 索引是区域编号，值代表这个区域方程在funcs方程中的位置
        end
        
        
        function cal_BC(obj, Ln, Rn, Tn, Bn)
            syms x y real
            
            BCfuncs_loc_num = 1;    % 用于记录当前边界函数在整体方程组中的位置
            
            
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
                            [~, eqT] = obj.genBB(top_idx, 1); % 这里取得的eqT是一个cell数组，其中的每一项是边界方程的一部分，例如c这部分，
                            for j = 1:length(eqT)
                                eqT{j} = subs(eqT{j}, y, obj.case_impl.yt);
                            end
                            % 如果是分段函数，传入的T_funcs将会是一个更长的一维cell数组
                            obj.case_impl.T_funcs = [obj.case_impl.T_funcs, eqT];     % 这样子就成为了一个一维数组，每次传入其中的一行
                            obj.BCfuncs_loc(:,top_idx) = [1; BCfuncs_loc_num];
                            BCfuncs_loc_num = BCfuncs_loc_num + 1;
                            
                            % 处理ES，如果这个边界是ES
                            if obj.case_impl.ES_regions(top_idx)
                                edge_impl = obj.impl.all_regions{top_idx}.impl;
                                ES = subs(edge_impl.B_x_P, y, obj.case_impl.yt);
                                obj.case_impl.T_ESfuncs = [obj.case_impl.T_ESfuncs; ES];
                            end
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
                            [~, eqB] = obj.genBB(bottom_idx, 1); % 这里取得的eqT是一个cell数组，其中的每一项是边界方程的一部分，例如c这部分，
                            for j = 1:length(eqB)
                                eqB{j} = subs(eqB{j}, y, obj.case_impl.yl);
                            end
                            % 如果是分段函数，传入的T_funcs将会是一个更长的一维cell数组
                            obj.case_impl.B_funcs = [obj.case_impl.B_funcs, eqB];     % 这样子就成为了一个一维数组，每次传入其中的一行
                            obj.BCfuncs_loc(:,bottom_idx) = [2; BCfuncs_loc_num];
                            BCfuncs_loc_num = BCfuncs_loc_num + 1;
                            
                            % 处理ES，如果这个边界是ES
                            if obj.case_impl.ES_regions(bottom_idx)
                                edge_impl = obj.impl.all_regions{bottom_idx}.impl;
                                ES = subs(edge_impl.B_x_P, y, obj.case_impl.yl);
                                obj.case_impl.B_ESfuncs = [obj.case_impl.B_ESfuncs; ES];
                            end
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
                        obj.BCfuncs_loc(:,left_idx) = [3;BCfuncs_loc_num];
                        BCfuncs_loc_num = BCfuncs_loc_num + 1;
                        
                        % 处理ES，如果这个边界是ES
                        if obj.case_impl.ES_regions(left_idx)
                            edge_impl = obj.impl.all_regions{left_idx}.impl;
                            ES = subs(edge_impl.A_y_P, x, obj.case_impl.xl);
                            ES = ES * obj.case_impl.lambda_n;
                            obj.case_impl.L_ESfuncs = [obj.case_impl.L_ESfuncs; ES];
                        end
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
                        obj.BCfuncs_loc(:,right_idx) = [4;BCfuncs_loc_num];
                        BCfuncs_loc_num = BCfuncs_loc_num + 1;
                        
                        % 处理ES，如果这个边界是ES
                        if obj.case_impl.ES_regions(right_idx)
                            edge_impl = obj.impl.all_regions{right_idx}.impl;
                            ES = subs(edge_impl.A_y_P, x, obj.case_impl.xr);
                            ES = ES * obj.case_impl.lambda_n;
                            obj.case_impl.R_ESfuncs = [obj.case_impl.R_ESfuncs; ES];
                        end
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
                        obj.BCfuncs_loc(:,top_idx) = [1;BCfuncs_loc_num];
                        BCfuncs_loc_num = BCfuncs_loc_num + 1;
                        
                        % 处理ES，如果这个边界是ES
                        if obj.case_impl.ES_regions(top_idx)
                            edge_impl = obj.impl.all_regions{top_idx}.impl;
                            ES = subs(edge_impl.A_y_P, y, obj.case_impl.yt);
                            ES = ES * obj.case_impl.beta_h;
                            obj.case_impl.T_ESfuncs = [obj.case_impl.T_ESfuncs; ES];
                        end
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
                        obj.BCfuncs_loc(:,bottom_idx) = [2;BCfuncs_loc_num];
                        BCfuncs_loc_num = BCfuncs_loc_num + 1;
                        
                        % 处理ES，如果这个边界是ES
                        if obj.case_impl.ES_regions(bottom_idx)
                            edge_impl = obj.impl.all_regions{bottom_idx}.impl;
                            ES = subs(edge_impl.A_y_P, y, obj.case_impl.yl);
                            ES = ES * obj.case_impl.beta_h;
                            obj.case_impl.B_ESfuncs = [obj.case_impl.B_ESfuncs; ES];
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
                        obj.BCfuncs_loc(:,left_idx) = [3;BCfuncs_loc_num];
                        BCfuncs_loc_num = BCfuncs_loc_num + 1;
                        
                        % 处理ES，如果这个边界是ES
                        if obj.case_impl.ES_regions(left_idx)
                            edge_impl = obj.impl.all_regions{left_idx}.impl;
                            ES = subs(edge_impl.A_y_P, x, obj.case_impl.xl);
                            ES = ES * obj.case_impl.lambda_n;
                            obj.case_impl.L_ESfuncs = [obj.case_impl.L_ESfuncs; ES];
                        end
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
                        obj.BCfuncs_loc(:,right_idx) = [4;BCfuncs_loc_num];
                        BCfuncs_loc_num = BCfuncs_loc_num + 1;
                        
                        % 处理ES，如果这个边界是ES
                        if obj.case_impl.ES_regions(right_idx)
                            edge_impl = obj.impl.all_regions{right_idx}.impl;
                            ES = subs(edge_impl.A_y_P, x, obj.case_impl.xr);
                            ES = ES * obj.case_impl.lambda_n;
                            obj.case_impl.R_ESfuncs = [obj.case_impl.R_ESfuncs; ES];
                        end
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
                            eqT{i} = eqT{i} * obj.case_impl.beta_h;  % 乘以 beta_h
                        end
                        obj.case_impl.T_funcs = eqT;
                        obj.BCfuncs_loc(:,top_idx) = [1;BCfuncs_loc_num];
                        BCfuncs_loc_num = BCfuncs_loc_num + 1;
                        
                        % 处理ES，如果这个边界是ES
                        if obj.case_impl.ES_regions(top_idx)
                            edge_impl = obj.impl.all_regions{top_idx}.impl;
                            ES = subs(edge_impl.A_y_P, y, obj.case_impl.yt);
                            ES = ES * obj.case_impl.beta_h;
                            obj.case_impl.T_ESfuncs = [obj.case_impl.T_ESfuncs; ES];
                        end
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
                        obj.BCfuncs_loc(:,bottom_idx) = [2;BCfuncs_loc_num];
                        BCfuncs_loc_num = BCfuncs_loc_num + 1;
                        
                        % 处理ES，如果这个边界是ES
                        if obj.case_impl.ES_regions(bottom_idx)
                            edge_impl = obj.impl.all_regions{bottom_idx}.impl;
                            ES = subs(edge_impl.A_y_P, y, obj.case_impl.yl);
                            ES = ES * obj.case_impl.beta_h;
                            obj.case_impl.B_ESfuncs = [obj.case_impl.B_ESfuncs; ES];
                        end
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
                        obj.BCfuncs_loc(:,left_idx) = [3;BCfuncs_loc_num];
                        BCfuncs_loc_num = BCfuncs_loc_num + 1;
                        
                        % 处理ES，如果这个边界是ES
                        if obj.case_impl.ES_regions(left_idx)
                            edge_impl = obj.impl.all_regions{left_idx}.impl;
                            ES = subs(edge_impl.B_y_P, x, obj.case_impl.xl);
                            obj.case_impl.L_ESfuncs = [obj.case_impl.L_ESfuncs; ES];
                        end
                    else
                        obj.case_impl.L_funcs = {};
                        obj.case_impl.L_coeffs = {};
                    end
                    
                    % 右边界
                    if Rn == 0 && ~isempty(obj.right)
                        right_idx = obj.right(1);
                        
                        [coeffs, eqR] = obj.genBB(right_idx, 2);
                        for i = 1:length(eqL)
                            eqR{i} = subs(eqR{i}, x, obj.case_impl.xr);
                        end
                        obj.case_impl.R_funcs = eqR;
                        obj.BCfuncs_loc(:,right_idx) = [4;BCfuncs_loc_num];
                        BCfuncs_loc_num = BCfuncs_loc_num + 1;
                        
                        % 处理ES，如果这个边界是ES
                        if obj.case_impl.ES_regions(right_idx)
                            edge_impl = obj.impl.all_regions{right_idx}.impl;
                            ES = subs(edge_impl.B_y_P, x, edge_impl.xr);
                            obj.case_impl.R_ESfuncs = [obj.case_impl.R_ESfuncs; ES];
                        end
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
                    eqF{1} = F;
                    coeffs{1} = edge_impl.c_0x;
                    
                    F = subs(edge_impl.A_zx_expr, {edge_impl.c_hx, edge_impl.d_hx, edge_impl.c_0x, edge_impl.d_0x}, {0, 0, 0, 1});    % 只留下c_0x
                    eqF{3} = F;
                    coeffs{3} = edge_impl.d_0x;
                end
                
                if edge_region.Bn == 0  % 相邻区域的这个边界上有c_hx
                    % 如果有c_0x，说明这个需要考虑分段函数
                    % 分段函数直接放在最后两个位置，分别是c_0x和d_0x的表达式
                    F1 = subs(edge_impl.A_zx_expr, {edge_impl.c_hx, edge_impl.d_hx}, {1, 0});
                    if has_cd0x
                        F1 = subs(F1, {edge_impl.c_0x, edge_impl.d_0x}, {0, 0});    % 如果有c_0x,d_0x的话，置为0
                    end
                    eqF{2} = F1;
                    coeffs{2} = edge_impl.c_hx;
                end
                if edge_region.Tn == 0
                    F1 = subs(edge_impl.A_zx_expr, {edge_impl.c_hx, edge_impl.d_hx}, {0, 1});
                    if has_cd0x
                        F1 = subs(F1, {edge_impl.c_0x, edge_impl.d_0x}, {0, 0});    % 如果有c_0x,d_0x的话，置为0
                    end
                    eqF{4} = F1;
                    coeffs{4} = edge_impl.d_hx;
                end
                if edge_region.Ln == 0
                    F1 = subs(edge_impl.A_zy_expr, {edge_impl.e_ny, edge_impl.f_ny}, {1, 0});
                    eqF{5} = F1;
                    coeffs{5} = edge_impl.e_ny;
                end
                if edge_region.Rn == 0
                    F1 = subs(edge_impl.A_zy_expr, {edge_impl.e_ny, edge_impl.f_ny}, {0, 1});
                    eqF{6} = F1;
                    coeffs{6} = edge_impl.f_ny;
                end
                
                % 处理ES
                
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
                    eqF{1} = F;
                    coeffs{1} = edge_impl.c_0x;
                    
                    F = subs(F1, {edge_impl.c_hx, edge_impl.d_hx, edge_impl.c_0x, edge_impl.d_0x}, {0, 0, 0, 1});    % 只留下c_0x
                    eqF{3} = F;
                    coeffs{3} = edge_impl.d_0x;
                end
                if edge_region.Bn == 0  % 相邻区域的这个边界上有c_hx
                    % 如果有c_0x，说明这个需要考虑分段函数
                    % 分段函数直接放在最后两个位置，分别是c_0x和d_0x的表达式
                    F1 = obj.getBB_func(edge_impl, 1, cd_or_ef);
                    F1 = subs(F1, {edge_impl.c_hx, edge_impl.d_hx}, {1, 0});
                    if has_cd0x
                        F1 = subs(F1, {edge_impl.c_0x, edge_impl.d_0x}, {0, 0});    % 如果有c_0x,d_0x的话，置为0
                    end
                    eqF{2} = F1;
                    coeffs{2} = edge_impl.c_hx;
                end
                if edge_region.Tn == 0  % 相邻区域的这个边界上有d_hx
                    % 如果有c_0x，说明这个需要考虑分段函数
                    % 分段函数直接放在最后两个位置，分别是c_0x和d_0x的表达式
                    F1 = obj.getBB_func(edge_impl, 1, cd_or_ef);
                    F1 = subs(F1, {edge_impl.c_hx, edge_impl.d_hx}, {0, 1});
                    if has_cd0x
                        F1 = subs(F1, {edge_impl.c_0x, edge_impl.d_0x}, {0, 0});    % 如果有c_0x,d_0x的话，置为0
                    end
                    eqF{4} = F1;
                    coeffs{4} = edge_impl.d_hx;
                end
                if edge_region.Ln == 0
                    F1 = obj.getBB_func(edge_impl, 2, cd_or_ef);
                    F = subs(F1, {edge_impl.e_ny, edge_impl.f_ny}, {1,0});
                    eqF{5} = F;
                    coeffs{5} = edge_impl.e_ny;
                end
                
                if edge_region.Rn == 0
                    F1 = obj.getBB_func(edge_impl, 2, cd_or_ef);
                    F = subs(F1, {edge_impl.e_ny, edge_impl.f_ny}, {0,1});
                    eqF{6} = F;
                    coeffs{6} = edge_impl.f_ny;
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
            if( x_or_y == 1)    % 上下边界为x，左右为y
                if cd_or_ef == 1
                    F1 = edge_impl.B_x_x;% 只含c d项
                else
                    F1 = edge_impl.B_y_x;   % 只含c d项
                end
            else
                if cd_or_ef == 1
                    F1 = edge_impl.B_x_y;
                else
                    F1 = edge_impl.B_y_y;
                end
            end
        end
        
        function cal_ES(obj, Ln, Rn, Tn, Bn)
            switch obj.bc_type
                case BC_TYPE.BBAA
                    
                case BC_TYPE.AAAA
                case BC_TYPE.AABB
            end
            
            
        end
    end
end
