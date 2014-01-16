function xWebNode = defineWebNodes(BLADE, WEB)

%% re-assign some structure variable names (for convenience)
NUM_SEC  = BLADE.NUM_SEC;
zSec     = BLADE.zSec;
chord    = BLADE.chord;
pitAxis  = BLADE.pitAxis;
NUM_WEBS = WEB.NUM_WEBS;
inbStn   = WEB.inbStn;
oubStn   = WEB.oubStn;
inbChLoc = WEB.inbChLoc;
oubChLoc = WEB.oubChLoc;

%%
xWebNode = NaN(NUM_SEC, NUM_WEBS);  % x/c location of the webs, use NaN as a placeholder where webs do not exist  
for n = 1:NUM_WEBS
    x  = zSec(inbStn(n):oubStn(n));
    pa = pitAxis(inbStn(n):oubStn(n));
    c  = chord(inbStn(n):oubStn(n));
    y1 = (inbChLoc(n) - pa(1))*c(1);
    y2 = (oubChLoc(n) - pa(end))*c(end);
    x1 = x(1);
    x2 = x(end);
    m  = (y2 - y1)/(x2 - x1);
    y  = m.*(x - x1) + y1;
    
    xWebNode(inbStn(n):oubStn(n),n) = y./c + pa;
end

for i = 1:NUM_SEC   
    if any( diff(xWebNode(i,:)) <= 0 )
        error(['ERROR: The webs are not allowed to cross, but they appear to cross at station ' num2str(i)]);
    end    
end

end % function defineWebNodes

