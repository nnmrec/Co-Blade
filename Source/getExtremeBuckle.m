function maxRb = getExtremeBuckle(BLADE, Buckle, OPT)

station = zeros(BLADE.NUM_SEC, 1);
root    = zeros(BLADE.NUM_SEC, 1);
spar    = zeros(BLADE.NUM_SEC, 1);
LEP     = zeros(BLADE.NUM_SEC, 1);
TEP     = zeros(BLADE.NUM_SEC, 1);
tip     = zeros(BLADE.NUM_SEC, 1);
webs    = zeros(BLADE.NUM_SEC, 1);
for i = 1:BLADE.NUM_SEC  
    station(i) = max([Buckle(i).Top; Buckle(i).Bot; Buckle(i).Web]);
    
    if i < OPT.INB_STN
        % blade root
        root(i) = max([Buckle(i).Top(1); Buckle(i).Bot(1)]); 
        
    elseif i >= OPT.INB_STN && i <= OPT.OUB_STN
        % TEP, spar cap, LEP, and webs
        LEP(i)  = max([Buckle(i).Top(1); Buckle(i).Bot(1)]);
        spar(i) = max([Buckle(i).Top(2); Buckle(i).Bot(2)]);
        TEP(i)  = max([Buckle(i).Top(3); Buckle(i).Bot(3)]);
        webs(i) = max([Buckle(i).Web]);
        
    elseif i > OPT.OUB_STN
        % blade tip
        tip(i) = max([Buckle(i).Top(1); Buckle(i).Bot(1)]);
        
    else
        error('[ERROR]: unrecognized number of panels')
    end
      
end
maxRb.root  = max(root);
maxRb.LEP   = max(LEP);
maxRb.spar  = max(spar);
maxRb.TEP   = max(TEP);
maxRb.tip   = max(tip);
maxRb.webs  = max(webs);
maxRb.blade = max(station(1:BLADE.NUM_SEC-OPT.OUB_STN)); % removes data near the blade tip which sometimes blows up (numerical blowup)
