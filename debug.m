Tx1 = 0.28;
Ty1 = 0.1;
Tx3 = 0.1;
Ty3 = 0.04;
beta3_h3 = pi / Tx3 ;
beta1_h1 = pi / Tx1;

syms x real;

val = (2 / Tx1) * coth(beta3_h3 * Ty3) * int(sin(beta3_h3*(x-x1))*sin(beta1_h1*(x-x1)), x1, Tx3);

disp(double(val));
