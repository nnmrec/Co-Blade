function [x_web y_web] = webCoords(Panel, xw_top, yw_top, xw_bot, yw_bot, nWebs, WEB_NODES)

nPanelsTop = Panel.Top.nPanels;
nPanelsBot = Panel.Bot.nPanels;

% compute the x-y coordinates of the inside edge of the top panels
x_t_ou = cell(nPanelsTop, 1);
y_t_ou = cell(nPanelsTop, 1);
x_t_in = cell(nPanelsTop, 1);
y_t_in = cell(nPanelsTop, 1);
for k = 1:nPanelsTop
    nPts   = numel(Panel.Top.x{k});
    x_t_ou{k} = Panel.Top.x{k}(1:nPts/2);
    y_t_ou{k} = Panel.Top.y{k}(1:nPts/2);
    x_t_in{k} = Panel.Top.x{k}(end:-1:nPts/2+1);
    y_t_in{k} = Panel.Top.y{k}(end:-1:nPts/2+1);    
end
x_top_ou = cat(1, x_t_ou{:});
y_top_ou = cat(1, y_t_ou{:});
x_top_in = cat(1, x_t_in{:});
y_top_in = cat(1, y_t_in{:});

% compute the x-y coordinates of the inside edge of the bottom panels
x_b_ou = cell(nPanelsBot, 1);
y_b_ou = cell(nPanelsBot, 1);
x_b_in = cell(nPanelsBot, 1);
y_b_in = cell(nPanelsBot, 1);
for k = 1:nPanelsBot
    nPts   = numel(Panel.Bot.x{k});
    x_b_ou{k} = Panel.Bot.x{k}(1:nPts/2);
    y_b_ou{k} = Panel.Bot.y{k}(1:nPts/2);
    x_b_in{k} = Panel.Bot.x{k}(end:-1:nPts/2+1);
    y_b_in{k} = Panel.Bot.y{k}(end:-1:nPts/2+1);    
end
x_bot_ou = cat(1, x_b_ou{:});
y_bot_ou = cat(1, y_b_ou{:});
x_bot_in = cat(1, x_b_in{:});
y_bot_in = cat(1, y_b_in{:});

% compute the x-y coordinates of the shear web midline
x_web = cell(nWebs, 1);
y_web = cell(nWebs, 1);
for i = 1:nWebs
            
    r_top = hypot(xw_top(i) - x_top_ou, yw_top(i) - y_top_ou);
    r_bot = hypot(xw_bot(i) - x_bot_ou, yw_bot(i) - y_bot_ou);
    [unused i_top] = min(r_top);
    [unused i_bot] = min(r_bot);
    
    x_web{i} = linspace(x_top_in(i_top), x_bot_in(i_bot), WEB_NODES)';
    y_web{i} = linspace(y_top_in(i_top), y_bot_in(i_bot), WEB_NODES)';
    
end


end % function webCoords