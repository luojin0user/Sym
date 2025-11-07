syms n_3 n_6 y real
f_ny_3 = int(-(n_3*sin((pi*n_3*(y - 100))/40)*sin((pi*n_6*(y - 100))/40))/(n_6*sinh((pi*n_6)/2)), y, 100, 140, 'Hold', true)/20
f_ny_3 = subs(f_ny_3, {n_3,n_6}, {1,1});
f_ny_3 = release(f_ny_3);
A = vpa(f_ny_3);