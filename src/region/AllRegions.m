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
            obj.regions{1} = Region(1, 1, 0, 28, 0, 10, obj.regions);
            obj.regions{1}.set_H_N_max(10, 12);
            obj.regions{1}.set_boundary(0, 0, 1, 0);
            K1 = obj.regions{1}.get_region_func();
            
            disp('K1(1):');
            pretty(K1{1});  % MATLAB 版本的 pprint
            
            % disp('BC1(1):');
            % pretty(BC1{1});
            
            
            %% 区域 3
            obj.regions{3} = Region(1, 3, 0, 10, 10, 14, obj.regions);
            obj.regions{3}.set_H_N_max(10, 12);
            obj.regions{3}.set_boundary(0, 1, 1, 1);
            obj.regions{3}.boundarys.bottom = 1;
            K3 = obj.regions{3}.get_region_func();
            
            obj.regions{3}.region_bc_func = obj.regions{3}.boundarys.cal_BC(0, 1, 1, 1);
            
            disp('K3(1):');
            pretty(K3{1});
            
            BC3 = obj.regions{3}.gen_region_BC_func();
            
            disp('BC3(1):');
            pretty(BC3{1});
            
            
            % 在所有的区域内部方程完成运算之后，再调用函数进行所有边界的运算
            
            %% 符号代入计算
            eq_main = lhs(BC3{1}) - rhs(BC3{1});
            eq_sub1 = K3{1};
            disp(eq_main)
            disp(eq_sub1)
            
            
            % 将 A_z_3 代入 d_hx_3_eq
            eq_final = subs(eq_main, lhs(eq_sub1), rhs(eq_sub1));
            
            % 去掉 == ，直接得到左边 - 右边形式
            
            % eq_tmp = subs(eq_sub1, lhs(eq_sub2), rhs(eq_sub2));
            
            % eq_final = subs(eq_main, lhs(eq_sub1), rhs(eq_sub2));
            
            % eq_final = simplify(expand(eq_final), 'Steps', 200, 'IgnoreAnalyticConstraints', true);
            
            disp('代入后的符号方程 eqc:');
            pretty(eq_final);
        end
    end
end
