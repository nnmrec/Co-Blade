function f = finiteDiff(x)

N    = length(x);
f    = zeros(N,1);
f(1) = x(2) - x(1);      % forward  difference
f(N) = x(N) - x(N-1);    % backward difference
for i = 2:N-1    
    % average the forward and backward differences
    f(i) = (x(i+1) - x(i-1)) / 2;          
end

end % function finiteDiff