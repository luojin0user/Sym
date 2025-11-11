Tx1 = 0.28;
Ty1 = 0.1;
Tx3 = 0.1;
Ty3 = 0.04;
beta3_h3 = pi / Tx3 ;
beta1_h1 = pi / Tx1;
lambda_n3 = 2*pi / Ty3;

y2=0.1;
y3=0.14;

syms x y real;

val = -2/Ty3 * lambda_n3 * ( y3*int( sin(lambda_n3 * (y - y2)) , y, y2, y3) - ...
    int( y * sin(lambda_n3 * (y - y2)), y, y2, y3)  ...
    );

disp(double(val));
