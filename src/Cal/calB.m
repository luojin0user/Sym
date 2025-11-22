function [Bx_new, By_new, Bz_new] = calB(x0, y0, z0)
% 当运行main后，保存了obj文件以及ICs文件后，可以直接使用这个文件绘图
% 前提是激励不变，如果改变激励，例如改变线圈位置、电流大小等，需要重新运行main函数

ICs_xoz = load("./mat/xoz/ICs.mat", 'ICs').ICs;
obj_xoz = load("./mat/xoz/obj.mat", 'obj').obj;

ICs_yoz = load("./mat/yoz/ICs.mat", 'ICs').ICs;
obj_yoz = load("./mat/yoz/obj.mat", 'obj').obj;

coil = load("./mat/RCoil.mat", 'obj').obj;

factor_x = coil.fx;
factor_y = coil.fy;

[Bx, Bz_xoz] = obj_xoz.cal_Bx_By(ICs_xoz, x0', z0');
[By, Bz_yoz] = obj_yoz.cal_Bx_By(ICs_yoz, y0', z0');

fx = factor_x(x0, y0, z0);
fy = factor_y(x0, y0, z0);

Bx_new = Bx .* fx .* sqrt(2);
Bz_xoz_new = Bz_xoz .* fx .* sqrt(2);

By_new = By .* fy .* sqrt(2);
Bz_yoz_new = Bz_yoz .* fy .* sqrt(2);

Bz_new = Bz_xoz_new + Bz_yoz_new;
end

