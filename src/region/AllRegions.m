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
            obj.regions{1} = Region(1, 1, 0, 28, 0, 10, 0, 0, 1, 0, 10, 12, obj.regions);
            obj.regions{1}.set_boundary();
            K1 = obj.regions{1}.get_region_solution_func();
            
            disp('K1(1):');
            pretty(K1{1});  % MATLAB 版本的 pprint
            
            % disp('BC1(1):');
            % pretty(BC1{1});
            
            
            %% 区域 3
            obj.regions{3} = Region(1, 3, 0, 10, 10, 14, 0, 1, 1, 1, 10, 12, obj.regions);
            obj.regions{3}.set_boundary();
            K3 = obj.regions{3}.get_region_solution_func();
            
            BC3 = obj.regions{3}.gen_region_coefficient_func();
            
            disp('K3(1):');
            pretty(K3{1});
            
            disp('BC3(1):');
            pretty(BC3{1}{1});
            
            
            % 在所有的区域内部方程完成运算之后，再调用函数进行所有边界的运算
            
            %% 符号代入计算
            eq_main = BC3{1}{1};
            eq_sub1 = K3{1};
            % disp(eq_main)
            % disp(eq_sub1)
            
            
            % 将 A_z_3 代入 d_hx_3_eq
            syms x y real
            eq_final = subs(eq_main, lhs(eq_sub1), rhs(eq_sub1));
            eq_final = subs(eq_final, obj.regions{3}.impl.h, 10);
            % 去除 Hold，自动求积分
            eq_final = release(rhs(eq_final));
            
            eq_final_value = double(eq_final);
            
            % 提取 d_hx_1 的系数（就是常数 A）
            % 去掉 == ，直接得到左边 - 右边形式
            
            % eq_tmp = subs(eq_sub1, lhs(eq_sub2), rhs(eq_sub2));
            
            % eq_final = subs(eq_main, lhs(eq_sub1), rhs(eq_sub2));
            
            % eq_final = simplify(expand(eq_final), 'Steps', 200, 'IgnoreAnalyticConstraints', true);
            
            disp('代入后的符号方程 eqc:');
            pretty(eq_final);
            disp('代入后的符号方程 eqc 数值:');
            disp(eq_final_value);
        end
    end
end
