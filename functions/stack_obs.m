function [zz] = stack_obs(z,p, drop_last)
    if nargin == 2 || ~drop_last
        [r,c] = size(z);
        zz = zeros(r-p+1, c*p);
        j = 1;
        for i = 1:p
                zz(:,j:j+c-1) = z((p-i+1):r-i+1,:);
                j = j+c;
        end
    else
       [r,c] = size(z);
        zz = zeros(r-p, c*p);
        j = 1;
        for i = 1:p
                zz(:,j:j+c-1) = z((p-i+1):r-i,:);
                j = j+c;
        end     
    end
end