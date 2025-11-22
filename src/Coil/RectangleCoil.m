classdef RectangleCoil < handle
    properties
        % 一个线圈的所有属性
        % 把一个线圈从顶部往下看
        % 三维信息
        
        % x和y方向的外长度
        outer_len_x
        outer_len_y
        
        % x和y方向的内长度
        inner_len_x
        inner_len_y
        
        % 理论上而言，一半的线圈的宽度在x和y方向是相同的，
        % 也就是 outer_len_y - inner_len_y==outer_len_x - inner_len_x
        
        % x,y,z中心位置，z是底面位置
        x_loc
        y_loc
        z_loc
        
        mu0
        Nt              % 总匝数
        Nt_per_layer    %每一层的匝数
        layer           % 层数，
        wire_diameter   % 导线直径
        thick    % 厚度
        Ir   % 电流大小，电流默认是顺时针，流入为正
        
        % 用于计算的区域
        cal_xl
        cal_xr
        cal_yl
        cal_yr
        cal_zb
        cal_zt
        
        fx
        fy
        
        
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
        function obj = RectangleCoil(outer_len_x, outer_len_y, Nt_per_layer, layer, wire_diameter, Ir)
            obj.outer_len_x = outer_len_x;
            obj.outer_len_y = outer_len_y;
            obj.Nt_per_layer = Nt_per_layer;
            obj.layer = layer;
            obj.Ir = Ir;
            obj.wire_diameter = wire_diameter;
            
            obj.inner_len_x = outer_len_x - Nt_per_layer * wire_diameter * 2;
            obj.inner_len_y = outer_len_y - Nt_per_layer * wire_diameter * 2;
            
            obj.Nt = Nt_per_layer * layer;
            obj.thick = wire_diameter * layer;
        end
        
        % 生成所有区域，并进行计算
        function gen_all_regions(obj)
            plain = AllRegions();  % 创建对象
            
            plain.input_current_region("xoz", obj.x_cl1, obj.x_cl2, obj.z_b, obj.z_t, 1, obj.Jr_xl);
            plain.input_current_region("xoz", obj.x_cr1, obj.x_cr2, obj.z_b, obj.z_t, 1, obj.Jr_xr);
            
            plain.input_current_region("yoz", obj.y_cl1, obj.y_cl2, obj.z_b, obj.z_t, 1, obj.Jr_yl);
            plain.input_current_region("yoz", obj.y_cr1, obj.y_cr2, obj.z_b, obj.z_t, 1, obj.Jr_yr);
            
            plain.input_calculate_area(obj.cal_xl, obj.cal_xr, obj.cal_zb, obj.cal_zt, 1);
            plain.pre_process();
            plain.get_all_regions();  % 调用方法生成并显示所有区域
            
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
            
            obj.z_b = z;
            obj.z_t = z + obj.thick;
            
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
            
            obj.factor_x();
            obj.factor_y();
        end
        
        % 用于2D到3D的修正因子计算
        % 输入一个线圈
        % 返回一个匿名函数
        function factor_x(obj)
            mu_0 = obj.mu0;
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
            
            B_inf_x = sqrt((mu_0 .* Ip ./ g1_2pi)^2 + (mu_0 .* Ip ./ g2_2pi)^2 ...
                - 2 .* (mu_0 .* Ip ./ g1_2pi) .* (mu_0 .* Ip ./ g2_2pi) .* cos(pi - alpha1 - alpha2));
            
            B_f_1 = mu_0 .* Ip ./ (2 .* g1_2pi) .* (g4 ./ sqrt(g1^2 + g4^2) - g3 ./ sqrt(g1^2 + g3^2));
            B_f_2 = mu_0 .* Ip ./ (2 .* g2_2pi) .* (g4 ./ sqrt(g2^2 + g4^2) - g3 ./ sqrt(g2^2 + g3^2));
            
            B_f_x = sqrt(B_f_1^2 + B_f_2^2 ...
                - 2 .* B_f_1 .* B_f_2 .* cos(pi - alpha1 - alpha2));
            
            expr =  B_f_x ./ B_inf_x;
            obj.fx = matlabFunction(expr, 'Vars', {x, y, z});
        end
        
        
        function factor_y(obj)
            mu_0 = obj.mu0;
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
            
            B_inf_y = sqrt((mu_0 .* Ip ./ g1_2pi)^2 + (mu_0 .* Ip ./ g2_2pi)^2 ...
                - 2 .* (mu_0 .* Ip ./ g1_2pi) .* (mu_0 .* Ip ./ g2_2pi) .* cos(pi - alpha1 - alpha2));
            
            B_f_1 = mu_0 .* Ip ./ (2 .* g1_2pi) .* (g4 ./ sqrt(g1^2 + g4^2) - g3 ./ sqrt(g1^2 + g3^2));
            B_f_2 = mu_0 .* Ip ./ (2 .* g2_2pi) .* (g4 ./ sqrt(g2^2 + g4^2) - g3 ./ sqrt(g2^2 + g3^2));
            
            B_f_y = sqrt(B_f_1^2 + B_f_2^2 ...
                - 2 .* B_f_1 .* B_f_2 .* cos(pi - alpha1 - alpha2));
            
            expr = B_f_y ./ B_inf_y;
            
            obj.fy = matlabFunction(expr, 'Vars', {x, y, z});
        end
        
        
        function plot3D(obj, ax)
            % ============ 读取参数 ============
            ox = obj.outer_len_x;
            oy = obj.outer_len_y;
            ix = obj.inner_len_x;
            iy = obj.inner_len_y;
            t  = obj.thick;
            
            xc = obj.x_loc;
            yc = obj.y_loc;
            z0 = obj.z_loc;  % 底面坐标
            
            % ============ 外长方体的8个顶点 ============
            Xo = ox/2;  Yo = oy/2;
            Zo1 = z0;        % bottom
            Zo2 = z0 + t;    % top
            
            outer = [
                Xo,  Yo, Zo1;
                -Xo,  Yo, Zo1;
                -Xo, -Yo, Zo1;
                Xo, -Yo, Zo1;
                Xo,  Yo, Zo2;
                -Xo,  Yo, Zo2;
                -Xo, -Yo, Zo2;
                Xo, -Yo, Zo2
                ];
            
            % 平移到 (xc, yc)
            outer(:,1) = outer(:,1) + xc;
            outer(:,2) = outer(:,2) + yc;
            
            % ============ 内孔的8个顶点 ============
            Xi = ix/2;  Yi = iy/2;
            
            inner = [
                Xi,  Yi, Zo1;
                -Xi,  Yi, Zo1;
                -Xi, -Yi, Zo1;
                Xi, -Yi, Zo1;
                Xi,  Yi, Zo2;
                -Xi,  Yi, Zo2;
                -Xi, -Yi, Zo2;
                Xi, -Yi, Zo2
                ];
            
            inner(:,1) = inner(:,1) + xc;
            inner(:,2) = inner(:,2) + yc;
            
            % ============ 绘制外框 ============
            
            obj.drawBox(outer, [0 0.447 0.741], ax);    % 蓝色外壳
            obj.drawBox(inner, [1 1 1], ax);            % 白色内孔（中空）
            
            xlabel(ax, 'X (m)');
            ylabel(ax, 'Y (m)');
            zlabel(ax, 'Z (m)');
            title(ax, '3D Rectangle Coil');
            
            view(ax,3);
        end
        
        
        % ====== 辅助函数：绘制一个长方体（8个顶点）======
        function drawBox(obj, P, colorFace, ax)
            % P 为 8×3 矩阵，顺序如下：
            % 1~4: bottom face, 5~8: top face
            
            faces = [
                1 2 3 4;   % bottom
                5 6 7 8;   % top
                1 2 6 5;   % front
                2 3 7 6;   % left
                3 4 8 7;   % back
                4 1 5 8    % right
                ];
            
            patch(ax, 'Vertices',P,'Faces',faces,...
                'FaceColor',colorFace,'FaceAlpha',0.4,...
                'EdgeColor','k','LineWidth',1);
        end
        
        
    end
end