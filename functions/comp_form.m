function [A] = comp_form(B)
[r,c] = size(B);
A = [B; eye(c-r), zeros(c-r, r)];
end