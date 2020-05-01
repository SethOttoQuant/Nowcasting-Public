% Testing SKF vs Ksmooth

A = Spec.A;
Q = Spec.Q;
R = Spec.R;
HJ = Spec.HJ;
Y = X';

Z_0 = zeros(sA,1);

xx = eye(sA^2) - kron(A,A);
vQ = reshape(Q, (sA)^2, 1);
V_0 = xx\vQ;
V_0 = reshape(V_0,sA,sA); 

[Zsmooth, Vsmooth, ~, loglik] = runKF(Y, A, HJ, Q, diag(R), Z_0, V_0);

S.loglik