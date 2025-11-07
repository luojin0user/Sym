integrand_func = @(h_1,n_3)log( integral(@(x)sin((h_1.*x.*pi)./2.8e+2).*sinh((n_3.*x.*pi)./4.0e+1),0.0,1.0e+2)./(sinh(n_3.*pi.*(5.0./2.0)).*1.4e+2) );

Q = log(integrand_func(1,91));
disp(Q);