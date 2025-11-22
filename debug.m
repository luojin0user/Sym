syms beta_h1 beta_h2 x

% 1. 积分的一般情况 (Case A: beta_h1 ~= beta_h2)
I_int_A = (1/2) * ( sin((beta_h1 - beta_h2)/10) / (beta_h1 - beta_h2) - ...
    sin((beta_h1 + beta_h2)/10) / (beta_h1 + beta_h2) );

% 2. 积分的特殊情况 (Case B: beta_h1 == beta_h2)
% 注意：当 beta_h1 == beta_h2 时，beta_h1 - beta_h2 = 0，
%       且 beta_h1 + beta_h2 = 2*beta_h1
I_int_B = (1/2) * ( 1/10 - sin(2*beta_h1/10) / (2*beta_h1) );

% 3. 使用 piecewise 合并表达式
% 语法: piecewise(condition_1, expression_1, condition_2, expression_2, ..., otherwise_expression)
I_int_piecewise = piecewise(beta_h1 == beta_h2, I_int_B, ...
    I_int_A); % 否则 (otherwise) 默认使用 I_int_A

% 4. 最终结果 (加入常数项)
C = -(25/7) * coth(19*beta_h2/500);
I_total_sym = C * I_int_piecewise;

% 简化（可选，但通常推荐）
% I_total_sym = (50*int(-(cosh((19*beta_h2)/500)*sin(beta_h1*x)*sin(beta_h2*x))/sinh((19*beta_h2)/500), x, 0, 1/10, 'Hold', true))/7;
I_total_sym = simplify(I_total_sym);

B_H1_data = rand(200);
B_H2_data = rand(200);
digits(32);
% 假设 B_H1_data 和 B_H2_data 是你的数据向量
tic
I_sym_vector = subs(I_total_sym, {beta_h1, beta_h2}, {B_H1_data, B_H2_data});
final_results = vpa(I_sym_vector);
toc
disp(size(final_results));
save("test.mat", 'final_results');