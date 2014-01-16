function stop = gradientSearchOutput(x, ...
                                     optimValues, ...
                                     state, ...
                                     SIM, ...
                                     OPT, ...
                                     BLADE, ...
                                     WEB, ...
                                     AF, ...
                                     OUT, ...
                                     z_oub, ...
                                     z_CP)
            
persistent fig1 

stop       = false;

if optimValues.fval == Inf
    % a severe error occured
    return
end

%% assign the elements of the design vector into meaningful variable names
x_current = x(:);

%%                                    
switch state
    case 'init'
        if OUT.PLOT_OPT_ITER
            % create the figures
            fig1 = figure('name', 'Current Best Point', ...
                          'color', 'white', ...
                          'units', 'normalized',...
                          'outerposition', [0.1 0.1 0.8 0.8]);
        end
        
    otherwise
    
end

if OPT.WRITE_X_ITER
    % the user has requested to write the design variables to a text file
    fmt = [ repmat('%18.16f  ', 1, numel(x_current)), '\r\n' ];
    fid = fopen([SIM.rootDir filesep SIM.case '_iterX.out'], 'a');
    fprintf(fid, fmt, x_current);       
    fclose(fid);
end 

if OUT.PLOT_OPT_ITER
    plotBestPoint(fig1, x_current, OPT, BLADE, WEB, AF, z_oub, z_CP);      
end
    
end % function getColor
