function f = fitnessFunction(SIM, ANLS, OPT, BLADE, WEB, StrProps, LaminaSS, Buckle, Disp, Modes)

persistent init_fitness

%% re-assign some structure variable names (for convenience)
NUM_SEC      = BLADE.NUM_SEC;
zSec         = BLADE.zSec;
ROT_SPD      = BLADE.ROT_SPD;
nWebs        = WEB.nWebs;
N_MODES      = ANLS.N_MODES;
MIN_FREQ_SEP = OPT.MIN_FREQ_SEP;
MAX_TIP_D    = OPT.MAX_TIP_D;

%%
s11_fcT_top = cell(NUM_SEC, 1); 
s11_fcC_top = cell(NUM_SEC, 1); 
s22_fcT_top = cell(NUM_SEC, 1); 
s22_fcC_top = cell(NUM_SEC, 1); 
s12_fcS_top = cell(NUM_SEC, 1); 

s11_fcT_bot = cell(NUM_SEC, 1); 
s11_fcC_bot = cell(NUM_SEC, 1); 
s22_fcT_bot = cell(NUM_SEC, 1); 
s22_fcC_bot = cell(NUM_SEC, 1); 
s12_fcS_bot = cell(NUM_SEC, 1); 

s11_fcT_web = cell(NUM_SEC, 1); 
s11_fcC_web = cell(NUM_SEC, 1); 
s22_fcT_web = cell(NUM_SEC, 1); 
s22_fcC_web = cell(NUM_SEC, 1); 
s12_fcS_web = cell(NUM_SEC, 1); 

buckle_fc_top  = cell(NUM_SEC, 1); 
buckle_fc_bot  = cell(NUM_SEC, 1); 
buckle_fc_web  = cell(NUM_SEC, 1); 

for i = 1:NUM_SEC
        
    % ---- Top ----
    s11_fcT_top{i}   = cell2mat( cat(1, LaminaSS(i).Top.s_11_fc_T) );
    s11_fcC_top{i}   = cell2mat( cat(1, LaminaSS(i).Top.s_11_fc_C) );
    s22_fcT_top{i}   = cell2mat( cat(1, LaminaSS(i).Top.s_22_fc_T) );
    s22_fcC_top{i}   = cell2mat( cat(1, LaminaSS(i).Top.s_22_fc_C) );
    s12_fcS_top{i}   = cell2mat( cat(1, LaminaSS(i).Top.s_12_fc_S) );
    buckle_fc_top{i} = Buckle(i).Top;
    
    % ---- Bottom ----
    s11_fcT_bot{i}   = cell2mat( cat(1, LaminaSS(i).Bot.s_11_fc_T) );
    s11_fcC_bot{i}   = cell2mat( cat(1, LaminaSS(i).Bot.s_11_fc_C) );
    s22_fcT_bot{i}   = cell2mat( cat(1, LaminaSS(i).Bot.s_22_fc_T) );
    s22_fcC_bot{i}   = cell2mat( cat(1, LaminaSS(i).Bot.s_22_fc_C) );
    s12_fcS_bot{i}   = cell2mat( cat(1, LaminaSS(i).Bot.s_12_fc_S) );
    buckle_fc_bot{i} = Buckle(i).Bot;
    
    % ---- Web ----
    if nWebs(i) >= 1
        s11_fcT_web{i}   = cell2mat( cat(1, LaminaSS(i).Web.s_11_fc_T) );
        s11_fcC_web{i}   = cell2mat( cat(1, LaminaSS(i).Web.s_11_fc_C) );
        s22_fcT_web{i}   = cell2mat( cat(1, LaminaSS(i).Web.s_22_fc_T) );
        s22_fcC_web{i}   = cell2mat( cat(1, LaminaSS(i).Web.s_22_fc_C) );
        s12_fcS_web{i}   = cell2mat( cat(1, LaminaSS(i).Web.s_12_fc_S) );
        buckle_fc_web{i} = Buckle(i).Web;
    else
        s11_fcT_web{i}   = [];
        s11_fcC_web{i}   = [];
        s22_fcT_web{i}   = [];
        s22_fcC_web{i}   = [];
        s12_fcS_web{i}   = [];
        buckle_fc_web{i} = [];
    end
    
end
s11_fcT   = cell2mat( [s11_fcT_top; s11_fcT_bot; s11_fcT_web] );
s11_fcC   = cell2mat( [s11_fcC_top; s11_fcC_bot; s11_fcC_web] );
s22_fcT   = cell2mat( [s22_fcT_top; s22_fcT_bot; s22_fcT_web] );
s22_fcC   = cell2mat( [s22_fcC_top; s22_fcC_bot; s22_fcC_web] );
s12_fcS   = cell2mat( [s12_fcS_top; s12_fcS_bot; s12_fcS_web] );
buckle_fc = cell2mat( [buckle_fc_top; buckle_fc_bot; buckle_fc_web] );


% penalty factors for exceeding maximum stress
penalty_s11_fcT = max(1, max(s11_fcT));
penalty_s11_fcC = max(1, max(s11_fcC));
penalty_s22_fcT = max(1, max(s22_fcT));
penalty_s22_fcC = max(1, max(s22_fcC));
penalty_s12_fcS = max(1, max(s12_fcS));

% penalty factor for exceeding buckling stress
penalty_buckle  = max(1, max(buckle_fc));

% penalty factor for exceeding the maximum allowable tip deflection
penalty_deflect = max(1, Disp.tipDeflect/MAX_TIP_D);

% penalty factors for the blade rotation freq. being too close to the blade
% natural freqs.
if N_MODES >= 1
    bladeRotFreq = ROT_SPD / 60; % blade rotation frequency (Hz)
    penalty_freq = max(1, max( MIN_FREQ_SEP ./ abs( Modes.natFreq - bladeRotFreq )));
else
    penalty_freq = 1;
end

bladeMass = StrProps.bladeMass;

% fitness value    
f = prod( [bladeMass;
           penalty_s11_fcT^2;
           penalty_s11_fcC^2;
           penalty_s22_fcT^2;
           penalty_s22_fcC^2;
           penalty_s12_fcS^2;
           penalty_buckle^2;
           penalty_deflect^2; 
           penalty_freq^2] );
       
% is this the initial blade mass? THIS IS A HACK--better to have init_fitness as an input to this function
if isempty(init_fitness)
    init_fitness = StrProps.bladeMass;
end
f = f / init_fitness; 

% f                    
% bladeMass
% penalty_s11_fcT
% penalty_s11_fcC
% penalty_s22_fcT
% penalty_s22_fcC
% penalty_s12_fcS
% penalty_buckle
% penalty_deflect
% penalty_freq


% w = 8/init_bladeMass;
% w = 2/init_bladeMass;
% f = w*bladeMass + sum( [penalty_s11_fcT;
%                         penalty_s11_fcC;
%                         penalty_s22_fcT;
%                         penalty_s22_fcC;
%                         penalty_s12_fcS;
%                         penalty_buckle;
%                         penalty_deflect; 
%                         penalty_freq].^2 );

      
       
% f = w*bladeMass + sum( [penalty_s11_fcT;
%     penalty_s11_fcC;
%     penalty_s22_fcT;
%     penalty_s22_fcC;
%     penalty_s12_fcS;
%     penalty_buckle;
%     penalty_deflect; 
%     penalty_freq].^2 );
% init_bladeMass

% f = prod( [bladeMass;                         % this does not work very well
%            penalty_s11_fcT^2;
%            penalty_s11_fcC^2;
%            penalty_s22_fcT^2;
%            penalty_s22_fcC^2;
%            penalty_s12_fcS^2;
%            penalty_buckle^2;
%            penalty_deflect^2; 
%            penalty_freq^2] );

% f = bladeMass/8 * sum( [penalty_s11_fcT^2;    % this does not work very well
%                         penalty_s11_fcC^2;
%                         penalty_s22_fcT^2;
%                         penalty_s22_fcC^2;
%                         penalty_s12_fcS^2;
%                         penalty_buckle^2;
%                         penalty_deflect^2; 
%                         penalty_freq^2] );

% f = bladeMass/8 * sum( [penalty_s11_fcT;      % this does not work very well
%                         penalty_s11_fcC;
%                         penalty_s22_fcT;
%                         penalty_s22_fcC;
%                         penalty_s12_fcS;
%                         penalty_buckle;
%                         penalty_deflect; 
%                         penalty_freq] );

% f =               sum( [penalty_s11_fcT^2;    % just testing to see what this does
%                         penalty_s11_fcC^2;
%                         penalty_s22_fcT^2;
%                         penalty_s22_fcC^2;
%                         penalty_s12_fcS^2;
%                         penalty_buckle^2;
%                         penalty_deflect^2; 
%                         penalty_freq^2] );

if OPT.WRITE_F_ALL
    % the user has requested to write the fitness value and penalty factors to a text file
    fmt = [ '%1.16e  %1.16e  ' repmat('%18.16f  ', 1, 9) '\r\n' ];
    fid = fopen([SIM.outputDir filesep SIM.case '_allF.out'], 'a');
    fprintf(fid, fmt, f, ...
                      bladeMass, ...
                      penalty_s11_fcT, ...
                      penalty_s11_fcC, ...
                      penalty_s22_fcT, ...
                      penalty_s22_fcC, ...
                      penalty_s12_fcS, ...
                      penalty_buckle, ...
                      penalty_deflect, ...
                      penalty_freq, ...
                      Modes.natFreq(1));      
    fclose(fid);
end

end % function fitnessFunction