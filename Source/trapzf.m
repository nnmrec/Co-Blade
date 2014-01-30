function z = trapzf(x,y)
% perform numerical integration using the trapezoid rule
% this function is similar to the native MATLAB function trapz, but it has
% been modified for increased speed

m = numel(x);
if m == 1;
    z = 0;
    return
end

% make sure that x and y are column vectors
x = x(:);
y = y(:);

z = diff(x)' * (y(1:m-1) + y(2:m))/2;

end

