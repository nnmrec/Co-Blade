function x = planformLine(zSec, chord, pitAxis, inbChLoc, oubChLoc, INB_STN, OUB_STN)

z_inb = zSec(INB_STN);
z_oub = zSec(OUB_STN);

x_inb = (inbChLoc - pitAxis(INB_STN))*chord(INB_STN);
x_oub = (oubChLoc - pitAxis(OUB_STN))*chord(OUB_STN);

m = (x_oub - x_inb)/(z_oub - z_inb);

x = m.*(zSec - z_inb) + x_inb; % x-coordinates in the blade planform
  
end % function planformLine

