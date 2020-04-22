A = readmatrix("params/A.csv");
Q = readmatrix("params/Q.csv");
CC = readmatrix("params/HJ.csv");
R = readmatrix("params/R.csv");
Z_0 = readmatrix("params/Z_0.csv");
V_0 = readmatrix("params/V_0.csv");
C = readmatrix("params/H.csv");

writematrix(C, "params/H.csv");

writematrix(VVsmooth(:,:,50), "params/VVsmooth.csv")

jim = Vsmooth(:,:,203);
al = VVsmooth(:,:,203);