function [ind r] = findClosest(x, val)
% this function finds the index of the closest value within a vector
% x MUST be a sorted vector with unique values

N   = length(val);
ind = zeros(N,1);
r   = zeros(N,1);
for n = 1:N
    [r(n) ind(n)] = min( abs(x - val(n)) );
end

end % function findClosest
