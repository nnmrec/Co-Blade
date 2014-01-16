function [Iflap Iedge angFlap] = principalValues(Ix, Iy, Ixy, aeroTwst)

a       = (Iy + Ix)./2;
b       = sqrt( ((Iy - Ix)./2).^2 + Ixy.^2 );
Iflap   = a - b;
Iedge   = a + b;
% angFlap = (1/2) .* atan(2.*Ixy ./ (Iy - Ix)) .* (180/pi);

% these angles should all be the same value, except when the argument
% becomes too large or too negative the trig functions wrap around their
% domains in different ways
ang1 = (1/2) .* atan(2.*Ixy ./ (Iy - Ix)) .* (180/pi);
ang2 = (1/2) .* asin(Ixy./b) .* (180/pi);
ang3 = (1/2) .* acos((Iy - Ix)./(2.*b)) .* (180/pi);

% pick the angle closest to the aerodynamic twist angle as the flapwise angle
angles = [ang1 ang2 ang3];
[unused ii] = min(abs(angles - 10), [], 2);

N       = numel(aeroTwst);
angFlap = zeros(N, 1);
for n = 1:N
    angFlap(n) = angles(n,ii(n));
end

end % function principalValues
