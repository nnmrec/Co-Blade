function z = cumtrapzf(x,y)
% perform cumulative numerical integration using the trapezoid rule
% this function is similar to the native MATLAB function cumtrapz, but it has
% been modified for increased speed

m = numel(x);
if m == 1;
    z = 0;
    return
end

% make sure that x and y are column vectors
x = x(:);
y = y(:);

dt = diff(x,1,1)/2;
z  = [0; cumsum(dt .* (y(1:m-1,:) + y(2:m,:)),1)];

end


