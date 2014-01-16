function [X Y] = scaleAirfoilCoords(x, y, pitAxis, chord, aeroTwst)

%% Translate the x/c-coordinates so the pitch axis is at the origin
x = x - pitAxis;

%% Scale x/c - y/c coordinates by the chord length
x = x .* chord;
y = y .* chord;

%% rotate x - y airfoil coordinates by the aero twist angle
aeroTwst = aeroTwst * pi/180; % convert degrees to radians
r        = sqrt(x.^2 + y.^2);
theta    = atan2(y, x);       
X        = r .* cos(theta + aeroTwst);
Y        = r .* sin(theta + aeroTwst);

end % function scaleAirfoilCoords
