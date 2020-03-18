%--------------------------------------------------------------------------

% Helper matrix J for observation equation

function CJ = get_CJ(C,frq,is_diff,p)
[k,r] = size(C);
lags = frq;
lags(is_diff,:) = arrayfun(@(x)(2*x-1),frq(is_diff,:));
pp = max([lags;p]);
CJ = zeros(k, r*pp);
CJ(:,1:r) = C; %for high frequency data
idx = find(frq>1); %for low frequency data
for i = 1:size(idx,1)
    j = idx(i);
    if is_diff(j)
        CJ(j,1:r*(2*frq(j)-1)) = C(j,:)*kron([1:frq(j),(frq(j)-1):-1:1]/frq(j), eye(r));
    else
        CJ(j,1:r*frq(j)) = C(j,:)*kron(ones(1,frq(j))/frq(j),eye(r));
    end
end 
return
end