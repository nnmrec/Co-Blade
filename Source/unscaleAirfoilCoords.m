function [X Y] = unscaleAirfoilCoords(x, y, pitAxis, chord, aeroTwst)

%% rotate x - y airfoil coordinates by the aero twist angle
aeroTwst = aeroTwst * pi/180; % convert degrees to radians
r        = sqrt(x.^2 + y.^2);
theta    = atan2(y, x);       
x        = r .* cos(theta - aeroTwst);
y        = r .* sin(theta - aeroTwst);

%% Scale x/c - y/c coordinates by the chord length
x = x ./ chord;
Y = y ./ chord;

%% Translate the x/c-coordinates so the pitch axis is at the origin
X = x + pitAxis;

%% For robustness, make sure the x_c coordinates are between 0 and 1
X(X<0) = 0;
X(X>1) = 1;

end % function scaleAirfoilCoords
