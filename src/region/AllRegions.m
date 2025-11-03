classdef AllRegions < handle
    properties
        regions  % cell array 存储所有区域
        
    end
    
    methods
        function obj = AllRegions()
            obj.regions = cell(1,7);
            
        end
        
        function get_all_regions(obj)
            %% 区域 1
            obj.regions{1} = Region(1, CaseType.BTAir, 0, 280, 0, 100, ...
                BC_TYPE.AAAA, [3,6,5,7,4],[0],[0],[0], 1, 1);
            K1 = obj.regions{1}.get_region_solution_func();
            
            disp('K1(1):');
            pretty(K1{1});  % MATLAB 版本的 pprint
            
            %% 区域 2
            obj.regions{2} = Region(2, CaseType.BTAir, 0, 280, 140, 240, ...
                BC_TYPE.BBAA, [0], [3,6,5,7,4],[0],[0], 1, 1);
            K1 = obj.regions{2}.get_region_solution_func();
            
            disp('K1(1):');
            pretty(K1{1});  % MATLAB 版本的 pprint
            
            %% 区域 3
            obj.regions{3} = Region(3, CaseType.NormalAir, 0, 100, 100, 140, ...
                BC_TYPE.AAAA, [2],[1],[0],[6], 1, 1);
            K3 = obj.regions{3}.get_region_solution_func();
            
            disp('K3(1):');
            pretty(K3{1});  % MATLAB 版本的 pprint
            
            %% 区域6
            obj.regions{6} = Region(6, CaseType.FerriteCurrent, 100, 280, 100, 140, ...
                BC_TYPE.AABB, [2],[1],[3],[5], 1, 1);
            K6 = obj.regions{6}.get_region_solution_func();
            disp('K6(1):');
            pretty(K6{1});  % MATLAB 版本的 pprint
            
            % 在所有的区域内部方程完成运算之后，再调用函数进行所有边界的运算
            [BC3, M3] = obj.regions{3}.gen_region_coefficient_func(obj.regions);
            
            
            
            disp('BC3(1):');
            pretty(BC3{1}{2});
            
            disp('BC3(4):');
            pretty(BC3{4}{1});
            
            
            %% 符号代入计算
            eq_main = BC3{4}{1};
            % eq_main_means = M3{4}{1};
            
            % 将 A_z_3 代入 d_hx_3_eq
            % eq_final = subs(eq_main, lhs(eq_sub1), rhs(eq_sub1));
            eq_final = subs(eq_main, obj.regions{3}.impl.h, 1);
            eq_final = subs(eq_final, obj.regions{3}.impl.n, 1);
            % 去除 Hold，自动求积分
            eq_final = release(rhs(eq_final));
            
            eq_final_value = -double(eq_final); % 系数都需要取负号
            
            disp('代入后的符号方程 eqc:');
            pretty(eq_final);
            disp('代入后的符号方程 eqc 数值:');
            disp(eq_final_value);
            %disp("对应的系数：");
            %disp(eq_main_means);
        end
    end
end
