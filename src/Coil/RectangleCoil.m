classdef RectangleCoil < handle
    properties
        % 一个线圈的所有属性
        % 把一个线圈从顶部往下看
        % 三维信息
        thick    % 厚度
        % x和y方向的外长度
        outer_len_x
        outer_len_y
        
        % x和y方向的内长度
        inner_len_x
        inner_len_y
        
        % 理论上而言，一半的线圈的宽度在x和y方向是相同的，
        % 也就是 outer_len_y - inner_len_y==outer_len_x - inner_len_x
        
        % x,y,z中心位置
        x_loc
        y_loc
        z_loc
        
        mu0
        Nt   % 匝数
        Ir   % 电流大小，电流默认是顺时针，流入为正

        % 用于计算的区域
        cal_xl
        cal_xr
        cal_yl
        cal_yr
        cal_zb
        cal_zt
        
        
        % 此后的属性是利于分区和计算的属性，无需输入
        % x，y方向的四个位置，x和y方向的左右
        x_cl1
        x_cl2
        x_cr1
        x_cr2
        
        m_xl    % (x_cl1 + x_cl2)/2
        m_xr
        
        y_cl1
        y_cl2
        y_cr1
        y_cr2
        
        m_yl
        m_yr
        
        % z方向的位置
        z_b % 下
        z_t % 上
        
        Jr_xl   % 电流密度
        Jr_xr
        Jr_yl
        Jr_yr
    end
    
    methods
        function obj = RectangleCoil(thick, outer_len_x, inner_len_x, outer_len_y, inner_len_y, Nt, Ir)
            obj.thick = thick;
            obj.outer_len_x = outer_len_x;
            obj.outer_len_y = outer_len_y;
            obj.inner_len_x = inner_len_x;
            obj.inner_len_y = inner_len_y;
            obj.Nt = Nt;
            obj.Ir = Ir;
        end

        % 生成所有区域，并进行计算
        function gen_all_regions(obj)
            % xoz方向 Bx Bz
            plain1 = AllRegions("xoz");  % 创建对象
            plain1.input_calculate_area(obj.cal_xl, obj.cal_xr, obj.cal_zb, obj.cal_zt, 1);
            plain1.input_current_region(obj.x_cl1, obj.x_cl2, obj.z_b, obj.z_t, 1, obj.Jr_xl);
            plain1.input_current_region(obj.x_cr1, obj.x_cr2, obj.z_b, obj.z_t, 1, obj.Jr_xr);
            plain1.pre_process();
            plain1.get_all_regions();  % 调用方法生成并显示所有区域
            
            % yoz方向 By Bz
            plain2 = AllRegions("yoz");  % 创建对象
            plain2.input_calculate_area(obj.cal_yl, obj.cal_yr, obj.cal_zb, obj.cal_zt, 1);
            plain2.input_current_region(obj.y_cl1, obj.y_cl2, obj.z_b, obj.z_t, 1, obj.Jr_yl);
            plain2.input_current_region(obj.y_cr1, obj.y_cr2, obj.z_b, obj.z_t, 1, obj.Jr_yr);
            plain2.pre_process();
            plain2.get_all_regions();  % 调用方法生成并显示所有区域

            % 保存当前对象
            save('./mat/RCoil.mat', 'obj');
        end
        
        % 设置方形线圈的位置，x，y是中心位置，xoy坐标系下
        function set_Rcoil_loc(obj, x, y, z)
            obj.x_loc = x;
            obj.y_loc = y;
            obj.z_loc = z;
            
            obj.x_cl1 = x - obj.outer_len_x ./ 2;
            obj.x_cl2 = x - obj.inner_len_x ./ 2;
            obj.x_cr1 = x + obj.inner_len_x ./ 2;
            obj.x_cr2 = x + obj.outer_len_x ./ 2;
            
            obj.y_cl1 = y - obj.outer_len_y ./ 2;
            obj.y_cl2 = y - obj.inner_len_y ./ 2;
            obj.y_cr1 = y + obj.inner_len_y ./ 2;
            obj.y_cr2 = y + obj.outer_len_y ./ 2;
            
            obj.m_xl = (obj.x_cl1 + obj.x_cl2) ./ 2;
            obj.m_xr = (obj.x_cr1 + obj.x_cr2) ./ 2;
            obj.m_yl = (obj.y_cl1 + obj.y_cl2) ./ 2;
            obj.m_yr = (obj.y_cr1 + obj.y_cr2) ./ 2;
            
            obj.z_b = z - obj.thick ./ 2;
            obj.z_t = z + obj.thick ./ 2;
            
            obj.Jr_xl = obj.Nt * obj.Ir / (obj.thick .* (obj.outer_len_x - obj.inner_len_x) ./ 2);
            obj.Jr_xr = -obj.Jr_xl;
            obj.Jr_yl = obj.Nt * obj.Ir / (obj.thick .* (obj.outer_len_y - obj.inner_len_y) ./ 2);
            obj.Jr_yr = -obj.Jr_yl;
        end

        function set_calculate_area(obj, xl, xr, yl, yr, zb, zt, mu0)
            % 输入电流区域的坐标
            % 输入的xl,xr,yb,yt分别是左侧x坐标，右侧x坐标，下侧y坐标，上侧y坐标，mu_r指的是这个区域的相对磁导率，一般为1，
            obj.cal_xl = xl;
            obj.cal_xr = xr;
            obj.cal_yl = yl;
            obj.cal_yr = yr;
            obj.cal_zb = zb;
            obj.cal_zt = zt;

            obj.mu0 = mu0;
        end
        
        % 用于2D到3D的修正因子计算
        % 输入一个线圈
        % 返回一个匿名函数
        function f = factor_x(obj)
            mu0 = obj.mu0;
            Ip = obj.Ir;
            Zc = obj.thick;
            
            syms x y z real
            g1 = sqrt((x - obj.m_xl)^2 + (z - Zc)^2);
            g2 = sqrt((x - obj.m_xr)^2 + (z - Zc)^2);
            g3 = obj.y_cl1 - y;
            g4 = obj.y_cr2 - y;
            
            g1_2pi = g1 .* 2 .* pi;
            g2_2pi = g2 .* 2 .* pi;
            
            alpha1 = acos((x - obj.m_xl) ./ g1);
            alpha2 = acos((obj.m_xr - x) ./ g2);
            
            B_inf_x = sqrt((mu0 .* Ip ./ g1_2pi)^2 + (mu0 .* Ip ./ g2_2pi)^2 ...
                - 2 .* (mu0 .* Ip ./ g1_2pi) .* (mu0 .* Ip ./ g2_2pi) .* cos(pi - alpha1 - alpha2));
            
            B_f_1 = mu0 .* Ip ./ (2 .* g1_2pi) .* (g4 ./ sqrt(g1^2 + g4^2) - g3 ./ sqrt(g1^2 + g3^2));
            B_f_2 = mu0 .* Ip ./ (2 .* g2_2pi) .* (g4 ./ sqrt(g2^2 + g4^2) - g3 ./ sqrt(g2^2 + g3^2));
            
            B_f_x = sqrt(B_f_1^2 + B_f_2^2 ...
                - 2 .* B_f_1 .* B_f_2 .* cos(pi - alpha1 - alpha2));
            
            expr =  B_f_x ./ B_inf_x;
            f = matlabFunction(expr, 'Vars', {x, y, z});
        end


        function f = factor_y(obj)
            mu0 = obj.mu0;
            Ip = obj.Ir;
            Zc = obj.thick;
        
            % ---- 注意：此处故意交换了 x <-> y ----
            syms x y z real
        
            % 3、4 导线沿 x 轴，因此距离公式中 roles of x,y are swapped
            g1 = sqrt((y - obj.m_yl)^2 + (z - Zc)^2);
            g2 = sqrt((y - obj.m_yr)^2 + (z - Zc)^2);
            g3 = obj.x_cl1 - x;
            g4 = obj.x_cr2 - x;
        
            g1_2pi = g1 .* 2 .* pi;
            g2_2pi = g2 .* 2 .* pi;
        
            alpha1 = acos((y - obj.m_yl) ./ g1);
            alpha2 = acos((obj.m_yr - y) ./ g2);
        
            B_inf_y = sqrt((mu0 .* Ip ./ g1_2pi)^2 + (mu0 .* Ip ./ g2_2pi)^2 ...
                - 2 .* (mu0 .* Ip ./ g1_2pi) .* (mu0 .* Ip ./ g2_2pi) .* cos(pi - alpha1 - alpha2));
        
            B_f_1 = mu0 .* Ip ./ (2 .* g1_2pi) .* (g4 ./ sqrt(g1^2 + g4^2) - g3 ./ sqrt(g1^2 + g3^2));
            B_f_2 = mu0 .* Ip ./ (2 .* g2_2pi) .* (g4 ./ sqrt(g2^2 + g4^2) - g3 ./ sqrt(g2^2 + g3^2));
        
            B_f_y = sqrt(B_f_1^2 + B_f_2^2 ...
                - 2 .* B_f_1 .* B_f_2 .* cos(pi - alpha1 - alpha2));
        
            expr = B_f_y ./ B_inf_y;
        
            f = matlabFunction(expr, 'Vars', {x, y, z});
        end

    end
end