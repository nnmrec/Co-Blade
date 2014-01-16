function [stop options optchanged] = patternSearchOutput(optimvalues, ...
                                                         options, ...
                                                         flag, ...
                                                         SIM, ...
                                                         OPT, ...
                                                         BLADE, ...
                                                         WEB, ...
                                                         AF, ...
                                                         OUT, ...
                                                         z_oub, ...
                                                         z_CP)

% OPTIMVALUES is a structure containing information about the state of the optimization:
%            x: current point X 
%    iteration: iteration number
%         fval: function value 
%     meshsize: current mesh size 
%    funccount: number of function evaluations
%       method: method used in last iteration 
%       TolFun: tolerance on fval
%         TolX: tolerance on X
%
%   OPTIONS: Options structure used by PATTERNSEARCH.
%
%   FLAG: Current state in which OutPutFcn is called. Possible values are:
%         init: initialization state 
%         iter: iteration state
%    interrupt: intermediate state
%         done: final state
% 		
%   STOP: A boolean to stop the algorithm.
%
%   OPTCHANGED: A boolean indicating if the options have changed.
%
%	See also PATTERNSEARCH, GA, PSOPTIMSET, SEARCHFCNTEMPLATE



%%
persistent fig1 

stop       = false;
optchanged = false;

if optimvalues.fval == Inf
    % a severe error occured
    return
end

x_current = optimvalues.x(:);

%%                                    
switch flag
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

end % function patternSearchOutput


