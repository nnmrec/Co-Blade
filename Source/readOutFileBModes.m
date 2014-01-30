function Modes = readOutFileBModes(SIM, ANLS)

fid = fopen([SIM.outputDir filesep SIM.case '_BModes.out'],'r');
if fid == -1
    error(['ERROR: Could not locate and open file ' SIM.outputDir filesep SIM.case '_BModes.out']);
end

for n = 1:6
    line = fgetl(fid);
    if line == -1
        % note: occasionally BModes will crash for unknown reasons and
        % produce a blank output file. In this situtation fgetl() reads a
        % "end-of-file" character.  For now, just detect this error and
        % produce some placeholder values.
        fprintf(1, 'WARNING: BModes failed to execute properly. Continuing anyways...\n');  % FUCK
        Modes.natFreq    = Inf;
        Modes.span_loc   = [];
        Modes.flap_disp  = [];
        Modes.flap_slope = [];
        Modes.lag_disp   = [];
        Modes.lag_slope  = [];
        Modes.twist      = [];
        fclose(fid);
        return
    end
end
   
natFreq    = zeros(ANLS.N_MODES, 1);
span_loc   =  cell(ANLS.N_MODES, 1);
flap_disp  =  cell(ANLS.N_MODES, 1);
flap_slope =  cell(ANLS.N_MODES, 1);
lag_disp   =  cell(ANLS.N_MODES, 1);
lag_slope  =  cell(ANLS.N_MODES, 1);
twist      =  cell(ANLS.N_MODES, 1);
for i = 1:ANLS.N_MODES
    
    for n = 1:2
        fgetl(fid);
    end
    
    % read the line containing the natural frequency
    freqLine   = fgetl(fid);
    pat        = '[^\[\]Hz\s]+';
    x          = regexpi(freqLine, pat, 'match');
    natFreq(i) = str2double( x{7} );

    for n = 1:3
        fgetl(fid);
    end

    % read the element table
    array = readCellArray(fid, '%f %f %f %f %f %f', ANLS.N_ELEMS + 1); 
    span_loc{i}   = array{1};	
    flap_disp{i}  = array{2};	
    flap_slope{i} = array{3};	
    lag_disp{i}   = array{4};	
    lag_slope{i}  = array{5};	
    twist{i}      = array{6};
    
end

fclose(fid);

%% Collect the output
Modes.natFreq    = natFreq;
Modes.span_loc   = span_loc;
Modes.flap_disp  = flap_disp;
Modes.flap_slope = flap_slope;
Modes.lag_disp   = lag_disp;
Modes.lag_slope  = lag_slope;
Modes.twist      = twist;

end % function readOutFileBModes
