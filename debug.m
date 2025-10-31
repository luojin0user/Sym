function eq_A_z = gen_A_z_1()
syms x y h H_max d_hx_1
A_z_1(x,y) = (28*d_hx_1*symsum((sin((pi*h*x)/28)*sinh((pi*h*y)/28)) ...
    /(h*sinh((5*pi*h)/14)), h, 1, H_max))/pi;
A_z = symfun(A_z_1, [x,y]);
eq_A_z = symfun(sym(['A_z']), [x y]) == A_z;
end


function A_z_3 = gen_A_z_3(A_z_1)
syms x y
A_z_3(x,y) = A_z_1(x,y);
end


function eq = gen_d_hx_3_eq(A_z_3)
syms x y h d_hx_3
d_hx_3(x,y) = int(A_z_3*sin((pi*h*x)/10), x, 0, 10, 'Hold', true)/5;
eq = symfun(sym(['d_hx']), [x,y]) == d_hx_3;
end


function out = expand_d_hx_3(eq, A_z_3, A_z_1)

% 1) 把等式转成 lhs-rhs = 0 形式（不依赖 ==）
expr = lhs(eq) - rhs(eq);

% 2) 代入 A_z_3 → A_z_1
expr = subs(expr, A_z_3, A_z_1);

% 3) 展开，但不执行积分求值
out = expand(expr);
end


A_z_1 = gen_A_z_1();
A_z_3 = gen_A_z_3(A_z_1);
d_hx_3_eq = gen_d_hx_3_eq(A_z_3);

disp(A_z_1)
disp(A_z_3)
disp(d_hx_3_eq)

final_expr = expand_d_hx_3(d_hx_3_eq, A_z_3, A_z_1);
pretty(final_expr);
