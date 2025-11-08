t1 = load("BC.mat");
t2 = load("ES.mat");

BC = t1.BC;
ES = t2.ES;
BCx = sparse(BC);

IC = lsqr(BCx, ES, 1e-6, 100);


save("IC.mat",'IC');