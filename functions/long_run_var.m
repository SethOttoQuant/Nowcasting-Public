function V_0 = long_run_var(A,Q)
sA = size(A,1);
xx = eye(sA^2) - kron(A,A);
vQ = reshape(Q, (sA)^2, 1);
V_0 = xx\vQ;
V_0 = reshape(V_0,sA,sA); 
V_0 = (V_0 + V_0')/2; %reduce rounding error
end

