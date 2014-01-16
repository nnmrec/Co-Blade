function writeInpFileLaminate(SIM, BLADE, WEB, SECNODES, LamData)

%% re-assign some structure variable names (for convenience)
NUM_SEC  = BLADE.NUM_SEC;
NUM_WEBS = WEB.NUM_WEBS;
nWebs    = WEB.nWebs;
nSegsTop = SECNODES.nSegsTop;
nSegsBot = SECNODES.nSegsBot;
xNodeTop = SECNODES.xNodeTop;
xNodeBot = SECNODES.xNodeBot;

%%
date = datestr(now, 'mmmm dd, yyyy HH:MM AM');

for i = 1:NUM_SEC
    
    fid = fopen([SIM.laminateDir filesep SIM.case '_OPTIMAL_' num2str(i) '.lam'], 'w');
    
    fprintf(fid, 'Generated on %s by Co-Blade v%s \r\n', date, SIM.version); 
    fprintf(fid, 'This line is for user comments\r\n');
    
    fprintf(fid, '*************************** TOP SURFACE ****************************\r\n');
    fprintf(fid, '%12.0f           N_scts(1):  no of sectors on top surface\r\n', nSegsTop(i));
    fprintf(fid, '\r\n');
    fprintf(fid, 'normalized chord location of nodes defining airfoil sectors boundaries (xsec_node)\r\n');
    fmt = [repmat('%12.8f ', 1, numel(xNodeTop{i}) ), '\r\n'];
    fprintf(fid, fmt, xNodeTop{i});
    for n = 1:nSegsTop(i)
        fprintf(fid, '..................................................................\r\n');
        fprintf(fid, 'Sect_num     N_laminas\r\n');
        fprintf(fid, '%6.0f    %6.0f \r\n', n, LamData(i).Top.nLam(n));
        fprintf(fid, '\r\n');
        fprintf(fid, 'lamina    num of  thickness   fibers_direction  composite_material ID\r\n');
        fprintf(fid, 'number    plies   of ply (m)       (deg)               (-)\r\n');
        fprintf(fid, 'lam_num  N_plies    Tply         Tht_lam            Mat_id\r\n');
        for m = 1:LamData(i).Top.nLam(n)
            fprintf(fid, '%7.0f  %7.0f  %10.8f  %10.1f            %4.0f  %s\r\n', ...
                         LamData(i).Top.lamNum{n}(m), LamData(i).Top.nPlies{n}(m), ...
                         LamData(i).Top.tPly{n}(m),   LamData(i).Top.fibAng{n}(m), ...
                         LamData(i).Top.matID{n}(m),  LamData(i).Top.matName{n}{m});
        end
    end
    fprintf(fid, '..................................................................\r\n');
    
    fprintf(fid, '\r\n');
    
    fprintf(fid, '*************************** BOTTOM SURFACE *************************\r\n');
    fprintf(fid, '%12.0f           N_scts(2):  no of sectors on top surface\r\n', nSegsBot(i));
    fprintf(fid, '\r\n');
    fprintf(fid, 'normalized chord location of nodes defining airfoil sectors boundaries (xsec_node)\r\n');
    fmt = [repmat('%12.8f ', 1, numel(xNodeBot{i}) ), '\r\n'];
    fprintf(fid, fmt, xNodeBot{i});
    for n = 1:nSegsBot(i)
        fprintf(fid, '..................................................................\r\n');
        fprintf(fid, 'Sect_num     N_laminas\r\n');
        fprintf(fid, '%6.0f    %6.0f \r\n', n, LamData(i).Bot.nLam(n));
        fprintf(fid, '\r\n');
        fprintf(fid, 'lamina    num of  thickness   fibers_direction  composite_material ID\r\n');
        fprintf(fid, 'number    plies   of ply (m)       (deg)               (-)\r\n');
        fprintf(fid, 'lam_num  N_plies    Tply         Tht_lam            Mat_id\r\n');
        for m = 1:LamData(i).Bot.nLam(n)
            fprintf(fid, '%7.0f  %7.0f  %10.8f  %10.1f            %4.0f  %s\r\n', ...
                         LamData(i).Bot.lamNum{n}(m), LamData(i).Bot.nPlies{n}(m), ...
                         LamData(i).Bot.tPly{n}(m),   LamData(i).Bot.fibAng{n}(m), ...
                         LamData(i).Bot.matID{n}(m),  LamData(i).Bot.matName{n}{m});
        end
    end
    fprintf(fid, '..................................................................\r\n');
    
    fprintf(fid, '\r\n');
     
    if nWebs(i) > 1
        fprintf(fid, '*************************** WEBS **************************************\r\n');
        fprintf(fid, 'Laminae schedule for webs (input required only if webs exist at this section):\r\n');
        fprintf(fid, '\r\n');
        
        for n = 1:NUM_WEBS
            fprintf(fid, '..................................................................\r\n');
            fprintf(fid, 'web_num      N_laminas\r\n');
            fprintf(fid, '%6.0f    %6.0f \r\n', n, LamData(i).Web.nLam(n));
            fprintf(fid, '\r\n');
            fprintf(fid, 'lamina    num of  thickness   fibers_direction  composite_material ID\r\n');
            fprintf(fid, 'number    plies   of ply (m)       (deg)               (-)\r\n');
            fprintf(fid, 'lam_num  N_plies    Tply         Tht_lam            Mat_id\r\n');
            for m = 1:LamData(i).Web.nLam(n)
                fprintf(fid, '%7.0f  %7.0f  %10.8f  %10.1f            %4.0f  %s\r\n', ...
                         LamData(i).Web.lamNum{n}(m), LamData(i).Web.nPlies{n}(m), ...
                         LamData(i).Web.tPly{n}(m),   LamData(i).Web.fibAng{n}(m), ...
                         LamData(i).Web.matID{n}(m),  LamData(i).Web.matName{n}{m});
            end
        end
        
    end
    
    fclose(fid);
    
end

end % function writeInpFileLaminate

