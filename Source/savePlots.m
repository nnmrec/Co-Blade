function savePlots(fig, figTitle, SAVE_FIG_FMT)

if strncmp(SAVE_FIG_FMT,'-fig',4)
	saveas(fig, figTitle, 'fig')
else       
    fig_ops = textscan(SAVE_FIG_FMT, '%s','delimiter', ',');
    ops     = cell(length(fig_ops{:}), 1);
    for n = 1:length(fig_ops{:})
        ops{n} = strcat(', ''', fig_ops{:}(n), '''');
    end        
    expr = ['export_fig(fig, figTitle', char( strcat(ops{:}) ), ')'];
    eval(expr);
end

end % function savePlots

