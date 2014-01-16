function [WEB SECNODES LamData] = structOptimize(SIM, ANLS, OPT, ENV, BLADE, WEB, AF, MATS, Coord, OUT)

%% re-assign some structure variable names (for convenience)
INB_STN      = OPT.INB_STN;
TRAN_STN     = OPT.TRAN_STN;
OUB_STN      = OPT.OUB_STN;
NUM_CP       = OPT.NUM_CP;
NUM_SEC      = BLADE.NUM_SEC;
BLD_LENGTH   = BLADE.BLD_LENGTH;
HUB_RAD      = BLADE.HUB_RAD;
zSec         = BLADE.zSec;
chord        = BLADE.chord;
pitAxis      = BLADE.pitAxis;
NUM_WEBS     = WEB.NUM_WEBS;
NormAFcoords = AF.NormAFcoords;

%%
z_oub   = zSec(TRAN_STN:OUB_STN);                           % z-coordinates between the outboard station and the end of transition station                             
z_CP    = linspace(zSec(TRAN_STN), zSec(OUB_STN), NUM_CP)';	% z-coordinates of design variable control points
numVars = 3 + 5*NUM_CP + 4;                                 % number of design variables

nSegsTopBot = zeros(NUM_SEC, 1);
for i = 1:BLADE.NUM_SEC
    if i <= INB_STN || i > OUB_STN
        nSegsTopBot(i) = 1;
    else
        nSegsTopBot(i) = 3;
    end
end

%% prepare debugging output files, if requested
if OPT.WRITE_X_ITER || OPT.WRITE_X_ALL
    % the user has requested to write the design variables to a text file
    % prepare the header data for this file
    fmtX    = [ repmat('%18.18s  ', 1, numVars), '\r\n' ];
    headerX = {'w_cap_inb';
               'w_cap_oub';
               't_blade_root';
               cellstr( strcat('t_blade_skin', num2str((1:NUM_CP)')) );
               cellstr( strcat('t_cap_uni',    num2str((1:NUM_CP)')) );
               cellstr( strcat('t_cap_core',   num2str((1:NUM_CP)')) );
               cellstr( strcat('t_lep_core',   num2str((1:NUM_CP)')) );
               cellstr( strcat('t_tep_core',   num2str((1:NUM_CP)')) );
               cellstr( strcat('t_web_skin',   num2str((1:2)')) );
               cellstr( strcat('t_web_core',   num2str((1:2)')) )};
    unitsX  = {'(-)';
               '(-)';
               '(m)';
               cellstr( repmat('(m)', NUM_CP, 1) );
               cellstr( repmat('(m)', NUM_CP, 1) );
               cellstr( repmat('(m)', NUM_CP, 1) );
               cellstr( repmat('(m)', NUM_CP, 1) );
               cellstr( repmat('(m)', NUM_CP, 1) );
               '(m)';
               '(m)';
               '(m)';
               '(m)'};
    headerX = cat(1, headerX{:});       
    unitsX  = cat(1, unitsX{:}); 
    
    if OPT.WRITE_X_ITER
        fid = fopen([SIM.outputDir filesep SIM.case '_iterX.out'], 'w');
        fprintf(fid, fmtX, headerX{:});       
        fprintf(fid, fmtX, unitsX{:});
        fclose(fid);
    end
    
    if OPT.WRITE_X_ALL
        fid = fopen([SIM.outputDir filesep SIM.case '_allX.out'], 'w');
        fprintf(fid, fmtX, headerX{:});       
        fprintf(fid, fmtX, unitsX{:});
        fclose(fid);
    end
end

if OPT.WRITE_F_ALL
    % the user has requested to write the fitness values to a text file
    % prepare the header data of this file
    fmtF    = [ '%22.22s  %22.22s  ' repmat( '%18.18s  ', 1, 9 ), '\r\n' ];
    headerF = {'fitness_value', '(-)'; ...
               'blade_mass',    '(kg)'; ...
               'p1',            '(-)'; ...
               'p2',            '(-)'; ...  
               'p3',            '(-)'; ... 
               'p4',            '(-)'; ... 
               'p5',            '(-)'; ... 
               'p6',            '(-)'; ... 
               'p7',            '(-)'; ... 
               'p8',            '(-)'; ... 
               '1st_freq',      '(Hz)'}; 
    
    fid = fopen([SIM.outputDir filesep SIM.case '_allF.out'], 'w');
    fprintf(fid, fmtF, headerF{:,1});
    fprintf(fid, fmtF, headerF{:,2});
    fclose(fid);
end
       
%% assign the upper and lower bounds
LB = zeros(numVars, 1);
UB = zeros(numVars, 1);

LB(1) = 0.05;   % bounds for the inboard spar cap width (non-dimensional)
UB(1) = min([2*pitAxis(INB_STN), 2*(1-pitAxis(INB_STN)), 0.95]);

LB(2) = 0.05;   % bounds for the outboard spar cap width (non-dimensional)
UB(2) = min([2*pitAxis(OUB_STN), 2*(1-pitAxis(OUB_STN)), 0.95]);

LB(3) = 0.04 * sqrt((HUB_RAD + BLD_LENGTH)/40); % 1e-6; 	% bounds for the blade root material thickness
UB(3) = chord(1) ./ 4;

maxPanThick      = 0.10 * interp1(zSec, [NormAFcoords.maxPerThick]' .* chord, z_CP);
LB(4:3+NUM_CP*5) = 1e-6;   % bounds for the LEP, spar cap, and TEP lamina thicknesses
UB(4:3+NUM_CP*5) = repmat(maxPanThick, 5, 1);

maxWebThick        = 0.25 * interp1(zSec, chord ./ (NUM_WEBS+1), zSec([INB_STN, OUB_STN]));
LB(4+NUM_CP*5:end) = 1e-6; % bounds for the web lamina thicknesses
UB(4+NUM_CP*5:end) = repmat(maxWebThick, 2, 1);

%% assign the linear inequality constraints
numLinCon  = 1 + ...            % spar cap width monotonically decreasing
             1 + ...        	% root thickness must greater than blade skin thickness
             (NUM_CP-1)*5 + ... % lamina thicknesses are monotonically decreasing in the panels
             2;                 % lamina thicknesses are monotonically decreasing in the webs
             
A_t = zeros(NUM_CP-1, NUM_CP);
for i = 1:NUM_CP-1
    A_t(i, i)   = -1;
    A_t(i, i+1) = 1;
end
a       = cell(1, 5);
[a{:}]  = deal(A_t);
A_thick = blkdiag(a{:});

% linear inequality constraints of form Ax <= b
A = zeros(numLinCon, numVars);
b = zeros(numLinCon, 1);

A(1,[1 2]) = [-chord(INB_STN) chord(OUB_STN)];
A(2,[3 4]) = [-1 1];
A(3:2+size(A_thick,1), 4:3+size(A_thick,2)) = A_thick;      % ensures that LEP, spar cap, and TEP panel thickness is decreasing
A(numLinCon-1:end, numVars-3:end) = [-1 1 0 0; 0 0 -1 1];   % ensure that web thickness is descreaing

%% make the initial guess for the design variables
if OPT.READ_INITX
    % read initial guess from text file
    fid = fopen([SIM.optimDir filesep OPT.INITX_FILE], 'r');
    if fid == -1
        error(['ERROR: Could not locate and open file ' [SIM.optimDir filesep OPT.INITX_FILE]]);
    end
    w_cap_inb    = readScalar(fid);
    w_cap_oub    = readScalar(fid);
    t_blade_root = readScalar(fid);
    t_blade_skin = readVector(fid, NUM_CP);
    t_cap_uni    = readVector(fid, NUM_CP);
    t_cap_core   = readVector(fid, NUM_CP);
    t_lep_core   = readVector(fid, NUM_CP);
    t_tep_core   = readVector(fid, NUM_CP);
    t_web_skin   = readVector(fid, 2);
    t_web_core   = readVector(fid, 2);
    fclose(fid);
    % concatenate all the design variables into a single column vector
    xo = [w_cap_inb; ...        % scalar: (1 by 1)
          w_cap_oub; ...        % scalar: (1 by 1)       
          t_blade_root'; ...	% scalar: (1 by 1)
          t_blade_skin'; ...	% vector: (NUM_CP by 1)     
          t_cap_uni'; ...   	% vector: (NUM_CP by 1)
          t_cap_core'; ...  	% vector: (NUM_CP by 1)         
          t_lep_core'; ...  	% vector: (NUM_CP by 1)
          t_tep_core'; ...  	% vector: (NUM_CP by 1)
          t_web_skin'; ...  	% vector: (2 by 1)
          t_web_core'];     	% vector: (2 by 1)
        
else
    % create an initial guess automatically    
    % initial spar cap width (non-dimensionalized by the chord length)
    w_cap_inb   = 0.5;
    w_cap_oub   = 0.25;
 
    % get the airfoil maximum thickness at the max chord station and outboard station
    t_max_inb = NormAFcoords(TRAN_STN).maxPerThick * chord(TRAN_STN);
    t_max_oub = NormAFcoords(OUB_STN).maxPerThick  * chord(OUB_STN);
    
    % initial thickness of the root material
    t_blade_root = 0.08 * sqrt((HUB_RAD + BLD_LENGTH)/40);

    % initial thickness of the weave material
    t_blade_skin = 0.02 .* linspace(t_max_inb, t_max_oub, NUM_CP)';

    % initial spar cap lamina thicknesses
    t_cap_uni  = 0.04  .* linspace(t_max_inb, t_max_oub, NUM_CP)';
    t_cap_core = 0.001 .* linspace(t_max_inb, t_max_oub, NUM_CP)';

    % initial leading edge panel lamina thicknesses
    t_lep_core = 0.005 .* linspace(t_max_inb, t_max_oub, NUM_CP)';

    % initial trailing edge panel lamina thicknesses
    t_tep_core = t_lep_core;
       
    % initial total thickness of the web laminates                  
    tWeb_inb = 0.01 * 1/(NUM_WEBS + 1) * chord(INB_STN);
    tWeb_oub = tWeb_inb;
    
    % initial web lamina thicknesses
    t_web_skin = linspace(tWeb_inb, tWeb_oub, 2)';
    t_web_core = 0.25 * t_web_skin;

    % concatenate all the design variables into a single column vector
    xo = [w_cap_inb; ...    % scalar: (1 by 1)
          w_cap_oub; ...   	% scalar: (1 by 1)       
          t_blade_root; ...	% scalar: (1 by 1)
          t_blade_skin; ...	% vector: (NUM_CP by 1)     
          t_cap_uni; ...   	% vector: (NUM_CP by 1)
          t_cap_core; ...  	% vector: (NUM_CP by 1)         
          t_lep_core; ...  	% vector: (NUM_CP by 1)
          t_tep_core; ...  	% vector: (NUM_CP by 1)
          t_web_skin; ...  	% vector: (2 by 1)
          t_web_core];     	% vector: (2 by 1)
end

%% now test to make sure that the intial guess is correct size and satisfies all constraints    
if  size(xo, 1) == numVars && size(xo, 2) == 1
    % size if correct
else
    error('ERROR: The provided initial value is the incorrect size.');
end 
if all(xo >= LB)
    % lower bounds are satisfied
else   
    fprintf(1,'[xo, LB, xo >= LB]\r\n');
    disp([xo, LB, xo >= LB]);
    error('ERROR: The provided initial value violate the lower bounds.');
end
if all(xo <= UB)
    % upper bounds are satisfied
else
    fprintf(1,'[xo, UB, xo <= UB]\r\n');
    disp([xo, UB, xo <= UB]);
    error('ERROR: The provided initial value violates the upper bounds .');
end
if all(A*xo <= b)     
    % linear inequality constraints are satisfied
else
    disp(A*xo <= b);
    error('ERROR: The provided initial value violates the linear inequality constraints.');
end
%% function handle to the optimization fitness function
OPT_FUN = @(xo) structFitness(xo, SIM, ANLS, OPT, ENV, BLADE, WEB, AF, MATS, Coord, z_oub, z_CP, nSegsTopBot);
             
% get the initial fitness value (so we can later normalize based on this initial value)
init_fitness = OPT_FUN(xo);

%% run the optimization
switch OPT.OPT_METHOD
    case 'PS'
        %% initiate pattern search optimization
        % read the pattern search options from the text input file                                                      
        fid = fopen([SIM.optimDir filesep 'options-PatternSearch.inp'], 'r');
        if fid == -1
            error(['ERROR: Could not locate and open file ' [SIM.optimDir filesep 'options-PatternSearch.inp']]);
        end
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        TolMesh         = readScalar(fid);
        TolX            = readScalar(fid);
        TolFun          = readScalar(fid);
        TolBind         = readScalar(fid);
        MaxIter         = readScalar(fid);
        MaxFunEvals     = readScalar(fid);
        TimeLimit       = readScalar(fid);
        MeshContraction = readScalar(fid);
        MeshExpansion   = readScalar(fid);
        MeshAccelerator = readString(fid);
        MeshRotate      = readString(fid);
        InitialMeshSize = readScalar(fid);
        ScaleMesh       = readString(fid);
        MaxMeshSize     = readScalar(fid);
        PollMethod      = readString(fid);
        CompletePoll    = readString(fid);
        PollingOrder    = readString(fid);
        Display         = readString(fid);
        Cache           = readString(fid);
        CacheSize       = readScalar(fid);
        CacheTol        = readScalar(fid);
        fclose(fid);

        if OUT.PLOT_OPT_ITER
            plotFcns = {@psplotbestf, @psplotmeshsize ,@psplotfuncount, @psplotbestx};
        else
            % do not create any plots
            plotFcns = [];
        end

        % function handle for the custom output function
        psCustomOutput = @(optimvalues, options, flag) patternSearchOutput(optimvalues, ...
                                                                           options, ...
                                                                           flag, ...
                                                                           SIM, ...
                                                                           OPT, ...
                                                                           BLADE, ...
                                                                           WEB, ...
                                                                           AF, ...
                                                                           OUT, ...
                                                                           z_oub, ...
                                                                           z_CP);

        OptionsPS = psoptimset('TolMesh',         TolMesh, ...
                               'TolCon',          1e-6, ...
                               'TolX',            TolX, ...
                               'TolFun',          TolFun, ...
                               'TolBind',         TolBind, ...
                               'MaxIter',         MaxIter, ...
                               'MaxFunEvals',     MaxFunEvals, ...
                               'TimeLimit',       TimeLimit, ...
                               'MeshContraction', MeshContraction, ...
                               'MeshExpansion',   MeshExpansion, ...
                               'MeshAccelerator', MeshAccelerator, ...
                               'MeshRotate',      MeshRotate, ...
                               'InitialMeshSize', InitialMeshSize, ...
                               'ScaleMesh',       ScaleMesh, ...
                               'MaxMeshSize',     MaxMeshSize, ...
                               'InitialPenalty',  10, ...
                               'PenaltyFactor',   100, ...
                               'PollMethod',      PollMethod, ... 
                               'CompletePoll',    CompletePoll, ...
                               'PollingOrder',    PollingOrder, ...     
                               'SearchMethod',    [], ... 
                               'CompleteSearch',  'on', ...
                               'Display',         Display, ...
                               'OutputFcns',      psCustomOutput, ...
                               'PlotFcns',        plotFcns, ...
                               'PlotInterval',    1, ...
                               'Cache',           Cache, ...
                               'CacheSize',       CacheSize, ...
                               'CacheTol',        CacheTol, ...
                               'Vectorized',      'off', ...
                               'UseParallel',     'never');

        [xBest fval exitflag output] = patternsearch(OPT_FUN, xo, A, b, [], [], LB, UB, [], OptionsPS);         
       
    case 'GS'
        %% initiate gradient search optimization
        % read the gradient search options from the text input file                                                      
        fid = fopen([SIM.optimDir filesep 'options-GradientSearch.inp'], 'r');
        if fid == -1
            error(['ERROR: Could not locate and open file ' [SIM.optimDir filesep 'options-GradientSearch.inp']]);
        end
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        Algorithm           = readString(fid);
        Diagnostics         = readString(fid);
        Display             = readString(fid);
        FinDiffType         = readString(fid);
        ScaleProblem        = readString(fid);
        SubproblemAlgorithm = readString(fid);
        Hessian             = readString(fid);
        MaxIter             = readScalar(fid);
        TolX                = readScalar(fid);
        TolFun              = readScalar(fid);
        
        if OUT.PLOT_OPT_ITER
            plotFcns = {@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotconstrviolation,@optimplotstepsize,@optimplotfirstorderopt};
        else
            % do not create any plots
            plotFcns = [];
        end
        
        % function handle for the custom output function
        gsCustomOutput = @(x,optimValues,state) gradientSearchOutput(x, ...
                                                                     optimValues, ...
                                                                     state, ...
                                                                     SIM, ...
                                                                     OPT, ...
                                                                     BLADE, ...
                                                                     WEB, ...
                                                                     AF, ...
                                                                     OUT, ...
                                                                     z_oub, ...
                                                                     z_CP);
                                                                 
        OptionsGS = optimset('Algorithm',              Algorithm, ...
                             'AlwaysHonorConstraints', 'bounds', ...                       
                             'Diagnostics',            Diagnostics, ...
                             'Display',                Display, ...
                             'FinDiffType',            FinDiffType, ...
                             'ScaleProblem',           ScaleProblem, ...
                             'SubproblemAlgorithm ',   SubproblemAlgorithm, ... 
                             'Hessian',                Hessian, ... 
                             'UseParallel',            'never', ...
                             'MaxIter',                MaxIter, ...
                             'TolX',                   TolX, ...
                             'TolFun',                 TolFun, ...
                             'TypicalX',               xo, ...
                             'OutputFcn',              gsCustomOutput, ...   
                             'PlotFcns',               plotFcns); 

        [xBest fval exitflag output] = fmincon(OPT_FUN, xo, A, b, [], [], LB, UB, [], OptionsGS);
    
    case 'PSO'
        %% initiate particle swarm optimization
        % read the particle swarm options from the text input file                                                      
        fid = fopen([SIM.optimDir filesep 'options-ParticleSwarm.inp'], 'r');
        if fid == -1
            error(['ERROR: Could not locate and open file ' [SIM.optimDir filesep 'options-ParticleSwarm.inp']]);
        end
        fgetl(fid);
        fgetl(fid);
        fgetl(fid); 
        CognitiveAttraction = readScalar(fid); 
        ConstrBoundary      = readString(fid); 
        Display             = readString(fid); 
        FitnessLimit        = readScalar(fid); 
        Generations         = readScalar(fid);
        MaxFunEvals         = readScalar(fid);
        PopulationSize      = readScalar(fid); 
        SocialAttraction    = readScalar(fid); 
        StallGenLimit       = readScalar(fid); 
        TimeLimit           = readScalar(fid); 
        TolFun              = readScalar(fid); 
        VelocityLimit       = readScalar(fid); 
        fclose(fid);

        if OUT.PLOT_OPT_ITER
%             plotFcns = {@psoplotbestf, @psoplotscores};
            plotFcns = {@psoplotbestf};
        else
            do not create any plots
            plotFcns = [];
        end

        % function handle for the custom output function
        psoCustomOutput = @(optimvalues, options, flag) particleSwarmOutput(optimvalues, ...
                                                                            options, ...
                                                                            flag, ...
                                                                            SIM, ...
                                                                            OPT, ...
                                                                            BLADE, ...
                                                                            WEB, ...
                                                                            AF, ...
                                                                            OUT, ...
                                                                            z_oub, ...
                                                                            z_CP);

                                                                        
        % find a feasible initial population
        initPop = constrainedRandomVector(PopulationSize,LB,UB,A,b);
        initPop(:,1) = xo(:); % include our "good" initial guess
                                                                        
        OptionsPSO = psooptimset('AccelerationFcn',     @psoiterate, ... 
                                 'CognitiveAttraction', CognitiveAttraction, ...
                                 'ConstrBoundary',      ConstrBoundary, ...
                                 'Display',             Display, ...
                                 'DemoMode',            'off', ...
                                 'FitnessLimit',        FitnessLimit, ...
                                 'Generations',         Generations, ...
                                 'MaxFunEvals',         MaxFunEvals, ...
                                 'HybridFcn',           [], ...
                                 'InitialPopulation',   initPop', ...
                                 'InitialVelocities',   [], ...
                                 'PlotFcns',            plotFcns, ...
                                 'PlotInterval',        1, ...
                                 'PopInitRange',        [LB'; UB'], ...
                                 'PopulationSize',      PopulationSize, ...
                                 'PopulationType',      'doubleVector', ...
                                 'OutputFcns',          psoCustomOutput, ...
                                 'SocialAttraction',    SocialAttraction, ...
                                 'StallGenLimit',       StallGenLimit, ...
                                 'TimeLimit',           TimeLimit, ...
                                 'TolFun',              TolFun, ...
                                 'TolCon',              1e-6, ...
                                 'Vectorized',          'off', ...
                                 'VelocityLimit',       VelocityLimit);

        [xBest fval exitflag output population scores] = pso(OPT_FUN, numVars, A, b, [], [], LB', UB', [], OptionsPSO);         
        
        case 'GA'
        %% initiate genetic algorithm optimization
        
        if OUT.PLOT_OPT_ITER
            plotFcns = {@gaplotbestf, @gaplotbestindiv, @gaplotdistance, @gaplotmaxconstr, @gaplotscores};
        else
            % do not create any plots
            plotFcns = [];
        end

        % function handle for the custom output function
        gaCustomOutput = @(options, state, flag) geneticAlgorithmOutput(options, ...
                                                                        state, ...
                                                                        flag, ...
                                                                        SIM, ...
                                                                        OPT, ...
                                                                        BLADE, ...
                                                                        WEB, ...
                                                                        AF, ...
                                                                        OUT, ...
                                                                        z_oub, ...
                                                                        z_CP);
                                                                    
        % read the genetic algorithm options from the text input file  
        % NOTE: this is done in a different way for GA than other above
        % (in future versions, consider changing how input files are read to match this method)                             
        fid = fopen([SIM.optimDir filesep 'options-GeneticAlgorithm.inp'], 'r');
        if fid == -1
            error(['ERROR: Could not locate and open file ' [SIM.optimDir filesep 'options-GeneticAlgorithm.inp']]);
        end
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        readingInput = 'Reading';
        temp2 = blanks(0);
        while readingInput == 'Reading'
            temp1 = fgetl(fid);
            if temp1 == -1;
                break;
            end
            temp2 = [temp2 temp1];
        end;
        evalc(temp2);
        fclose(fid);
        
        % find a feasible initial population
        initPop = constrainedRandomVector(OptionsGA.PopulationSize,LB,UB,A,b);
        initPop(:,1) = xo(:); % include our "good" initial guess
        OptionsGA.InitialPopulation = initPop';
        
        % start the GA                           
        [xBest fval exitflag output population scores] = ga(OPT_FUN, numVars, A, b, [], [], LB', UB', [], OptionsGA);
              
end % switch OPT.OPT_METHOD

%% post-processing of optimization results
% optimal laminate data
[WEB SECNODES LamData] = defineLaminateData(xBest, SIM, OPT, BLADE, WEB, MATS, z_oub, z_CP, nSegsTopBot);

% write the optimal laminate data as text laminate input files
writeInpFileLaminate(SIM, BLADE, WEB, SECNODES, LamData);

% write the best values for the design variables
fid = fopen([SIM.outputDir filesep SIM.case '_bestX.out'], 'w');
fmt = [repmat('%1.16e  ', 1, numel(xBest)) '\r\n'];
fprintf(fid, fmt, xBest);
fclose(fid);

end % function structOptimize                                                 