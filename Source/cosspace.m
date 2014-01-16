function x = cosspace(x1, x2, N, flag)
%if flag = 'start', half cosing spacing, points are packed more dense at the beginning
%if flag = 'end',   half cosing spacing, points are packed more dense at the end
%if flag = 'both',  full cosing spacing, points are packed more dense at the beginning AND end

L = x2 - x1;

switch flag
    case 'start'
        t = linspace(pi,pi/2,N)';
        x = x1 + L.*(1 + cos(t));
    case 'end'
        t = linspace(pi/2,0,N)';
        x = x1 + L.*cos(t);
    case 'both'
        t = linspace(pi,0,N)';
        x = x1 + L.*(1 + cos(t))./2;
end

x(1)   = x1;    % fixes rounding errors from floating point
x(end) = x2;    % fixes rounding errors from floating point

end % function cosspace
