syms n_4 n_7 y 'real';
f_ny_7(n_4, n_7) = int((sin((pi.*n_4.*(y - 100))/40).*sin((pi.*n_7.*(y - 100))/40).*cosh((5.*pi.*n_4)/2))/sinh((5.*pi.*n_4)/2), y, 100, 140, 'Hold', true)/20;
a = simplifyFraction(f_ny_7)
b = simplify(f_ny_7)
% d = subs(a, {n_4, n_7}, {1,1})
