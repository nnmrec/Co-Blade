function maxPerBlade = getExtremeStress(matID, BLADE, Panel, LaminaSS)

maxPerTop.s_11_fc_T = nan(BLADE.NUM_SEC, 1);
maxPerTop.s_11_fc_C = nan(BLADE.NUM_SEC, 1);
maxPerTop.s_22_fc_T = nan(BLADE.NUM_SEC, 1);
maxPerTop.s_22_fc_C = nan(BLADE.NUM_SEC, 1);
maxPerTop.s_12_fc_S = nan(BLADE.NUM_SEC, 1);

maxPerBot.s_11_fc_T = nan(BLADE.NUM_SEC, 1);
maxPerBot.s_11_fc_C = nan(BLADE.NUM_SEC, 1);
maxPerBot.s_22_fc_T = nan(BLADE.NUM_SEC, 1);
maxPerBot.s_22_fc_C = nan(BLADE.NUM_SEC, 1);
maxPerBot.s_12_fc_S = nan(BLADE.NUM_SEC, 1);

maxPerWeb.s_11_fc_T = nan(BLADE.NUM_SEC, 1);
maxPerWeb.s_11_fc_C = nan(BLADE.NUM_SEC, 1);
maxPerWeb.s_22_fc_T = nan(BLADE.NUM_SEC, 1);
maxPerWeb.s_22_fc_C = nan(BLADE.NUM_SEC, 1);
maxPerWeb.s_12_fc_S = nan(BLADE.NUM_SEC, 1);

for i = 1:BLADE.NUM_SEC
    % top 
    perTopPanel.s_11_fc_T = nan(Panel(i).Top.nPanels, 1);
    perTopPanel.s_11_fc_C = nan(Panel(i).Top.nPanels, 1);
    perTopPanel.s_22_fc_T = nan(Panel(i).Top.nPanels, 1);
    perTopPanel.s_22_fc_C = nan(Panel(i).Top.nPanels, 1);
    perTopPanel.s_12_fc_S = nan(Panel(i).Top.nPanels, 1);
    for n = 1:Panel(i).Top.nPanels
        mask_matID = Panel(i).Top.matID{n} == matID;
        num_matID  = sum(mask_matID); 
        if num_matID >= 1
            s_11_fc_T = LaminaSS(i).Top(n).s_11_fc_T(mask_matID);
            s_11_fc_C = LaminaSS(i).Top(n).s_11_fc_C(mask_matID);
            s_22_fc_T = LaminaSS(i).Top(n).s_22_fc_T(mask_matID);
            s_22_fc_C = LaminaSS(i).Top(n).s_22_fc_C(mask_matID);
            s_12_fc_S = LaminaSS(i).Top(n).s_12_fc_S(mask_matID);          
            perTopPanel.s_11_fc_T(n) = max(max( [s_11_fc_T{:}] ));
            perTopPanel.s_11_fc_C(n) = max(max( [s_11_fc_C{:}] ));
            perTopPanel.s_22_fc_T(n) = max(max( [s_22_fc_T{:}] ));
            perTopPanel.s_22_fc_C(n) = max(max( [s_22_fc_C{:}] ));
            perTopPanel.s_12_fc_S(n) = max(max( [s_12_fc_S{:}] ));
        end       
    end
    maxPerTop.s_11_fc_T(i) = max(perTopPanel.s_11_fc_T);
    maxPerTop.s_11_fc_C(i) = max(perTopPanel.s_11_fc_C);
    maxPerTop.s_22_fc_T(i) = max(perTopPanel.s_22_fc_T);
    maxPerTop.s_22_fc_C(i) = max(perTopPanel.s_22_fc_C);
    maxPerTop.s_12_fc_S(i) = max(perTopPanel.s_12_fc_S);
   
    % bottom 
    perBotPanel.s_11_fc_T = nan(Panel(i).Bot.nPanels, 1);
    perBotPanel.s_11_fc_C = nan(Panel(i).Bot.nPanels, 1);
    perBotPanel.s_22_fc_T = nan(Panel(i).Bot.nPanels, 1);
    perBotPanel.s_22_fc_C = nan(Panel(i).Bot.nPanels, 1);
    perBotPanel.s_12_fc_S = nan(Panel(i).Bot.nPanels, 1);
    for n = 1:Panel(i).Bot.nPanels
        mask_matID = Panel(i).Bot.matID{n} == matID;
        num_matID  = sum(mask_matID); 
        if num_matID >= 1
            s_11_fc_T = LaminaSS(i).Bot(n).s_11_fc_T(mask_matID);
            s_11_fc_C = LaminaSS(i).Bot(n).s_11_fc_C(mask_matID);
            s_22_fc_T = LaminaSS(i).Bot(n).s_22_fc_T(mask_matID);
            s_22_fc_C = LaminaSS(i).Bot(n).s_22_fc_C(mask_matID);
            s_12_fc_S = LaminaSS(i).Bot(n).s_12_fc_S(mask_matID);          
            perBotPanel.s_11_fc_T(n) = max(max( [s_11_fc_T{:}] ));
            perBotPanel.s_11_fc_C(n) = max(max( [s_11_fc_C{:}] ));
            perBotPanel.s_22_fc_T(n) = max(max( [s_22_fc_T{:}] ));
            perBotPanel.s_22_fc_C(n) = max(max( [s_22_fc_C{:}] ));
            perBotPanel.s_12_fc_S(n) = max(max( [s_12_fc_S{:}] ));
        end       
    end
    maxPerBot.s_11_fc_T(i) = max(perBotPanel.s_11_fc_T);
    maxPerBot.s_11_fc_C(i) = max(perBotPanel.s_11_fc_C);
    maxPerBot.s_22_fc_T(i) = max(perBotPanel.s_22_fc_T);
    maxPerBot.s_22_fc_C(i) = max(perBotPanel.s_22_fc_C);
    maxPerBot.s_12_fc_S(i) = max(perBotPanel.s_12_fc_S);
    
    % webs 
    if Panel(i).Web.nPanels >= 1
        perWebPanel.s_11_fc_T = nan(Panel(i).Web.nPanels, 1);
        perWebPanel.s_11_fc_C = nan(Panel(i).Web.nPanels, 1);
        perWebPanel.s_22_fc_T = nan(Panel(i).Web.nPanels, 1);
        perWebPanel.s_22_fc_C = nan(Panel(i).Web.nPanels, 1);
        perWebPanel.s_12_fc_S = nan(Panel(i).Web.nPanels, 1);
        for n = 1:Panel(i).Web.nPanels
            mask_matID = Panel(i).Web.matID{n} == matID;
            num_matID  = sum(mask_matID); 
            if num_matID >= 1
                s_11_fc_T = LaminaSS(i).Web(n).s_11_fc_T(mask_matID);
                s_11_fc_C = LaminaSS(i).Web(n).s_11_fc_C(mask_matID);
                s_22_fc_T = LaminaSS(i).Web(n).s_22_fc_T(mask_matID);
                s_22_fc_C = LaminaSS(i).Web(n).s_22_fc_C(mask_matID);
                s_12_fc_S = LaminaSS(i).Web(n).s_12_fc_S(mask_matID);          
                perWebPanel.s_11_fc_T(n) = max(max( [s_11_fc_T{:}] ));
                perWebPanel.s_11_fc_C(n) = max(max( [s_11_fc_C{:}] ));
                perWebPanel.s_22_fc_T(n) = max(max( [s_22_fc_T{:}] ));
                perWebPanel.s_22_fc_C(n) = max(max( [s_22_fc_C{:}] ));
                perWebPanel.s_12_fc_S(n) = max(max( [s_12_fc_S{:}] ));
            end       
        end
        maxPerWeb.s_11_fc_T(i) = max(perWebPanel.s_11_fc_T);
        maxPerWeb.s_11_fc_C(i) = max(perWebPanel.s_11_fc_C);
        maxPerWeb.s_22_fc_T(i) = max(perWebPanel.s_22_fc_T);
        maxPerWeb.s_22_fc_C(i) = max(perWebPanel.s_22_fc_C);
        maxPerWeb.s_12_fc_S(i) = max(perWebPanel.s_12_fc_S);
    end
       
end

maxPerStation.s_11_fc_T = max([maxPerTop.s_11_fc_T maxPerBot.s_11_fc_T maxPerWeb.s_11_fc_T]');
maxPerStation.s_11_fc_C = max([maxPerTop.s_11_fc_C maxPerBot.s_11_fc_C maxPerWeb.s_11_fc_C]');
maxPerStation.s_22_fc_T = max([maxPerTop.s_22_fc_T maxPerBot.s_22_fc_T maxPerWeb.s_22_fc_T]');
maxPerStation.s_22_fc_C = max([maxPerTop.s_22_fc_C maxPerBot.s_22_fc_C maxPerWeb.s_22_fc_C]');
maxPerStation.s_12_fc_S = max([maxPerTop.s_12_fc_S maxPerBot.s_12_fc_S maxPerWeb.s_12_fc_S]');

maxPerBlade.s_11_fc_T   = max(maxPerStation.s_11_fc_T);
maxPerBlade.s_11_fc_C   = max(maxPerStation.s_11_fc_C);
maxPerBlade.s_22_fc_T   = max(maxPerStation.s_22_fc_T);
maxPerBlade.s_22_fc_C   = max(maxPerStation.s_22_fc_C);
maxPerBlade.s_12_fc_S   = max(maxPerStation.s_12_fc_S);




