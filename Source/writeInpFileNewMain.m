function writeInpFileNewMain(SIM, OPT, BLADE, WEB)

% create a copy of the main input file, but update any parameters that were changed/created by the optimization routine
fid1 = fopen(SIM.inputFile, 'r');
fid2 = fopen([SIM.rootDir filesep 'Inputs' filesep SIM.case '_OPTIMAL.inp'], 'w');
for n = 1:14
    line = fgetl(fid1);
    fprintf(fid2, '%s\r\n', line);
end
fgetl(fid1);
fprintf(fid2, 'false           OPTIMIZE:       Perform optimization of composite layup?\r\n');
line = fgetl(fid1);
fprintf(fid2, '%s\r\n', line);
fgetl(fid1);
fprintf(fid2, 'false           OPT_PITAXIS:    Optimize the pitch axis?\r\n');
for n = 1:32
    line = fgetl(fid1);
    fprintf(fid2, '%s\r\n', line);
end
for n = 1:BLADE.NUM_SEC
    fgetl(fid1);
    fprintf(fid2, '%5.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f    %s    %s\r\n', ...
                 BLADE.zFrac(n), BLADE.aeroTwst(n), BLADE.chord(n), BLADE.pitAxis(n), ...
                 BLADE.px_a(n),  BLADE.py_a(n),     BLADE.qz_a(n),  BLADE.afFile{n}, BLADE.strFile{n});
end
for n = 1:3
    line = fgetl(fid1);
    fprintf(fid2, '%s\r\n', line);
end
fgetl(fid1);
fgetl(fid1);
fprintf(fid2, 'webNum    inbStn    oubStn     inbChLoc     oubChLoc (This table of values is ignored if OPTIMIZE = true)\r\n');     
fprintf(fid2, '   (-)       (-)       (-)          (-)          (-)\r\n');
if OPT.OPTIMIZE
    for n = 1:WEB.NUM_WEBS
        fprintf(fid2, '%6.0f %9.0f %9.0f %12.8f %12.8f\r\n', n, OPT.INB_STN, OPT.OUB_STN, WEB.inbChLoc(n), WEB.oubChLoc(n));
    end
else
    for n = 1:WEB.NUM_WEBS
        line = fgetl(fid1);
        fprintf(fid2, '%s\r\n', line);
    end
end
line = [];
while ~strncmpi(line,'-----  Output Options', 21)
    line = fgetl(fid1);
end
fprintf(fid2,'-----  Output Options  ---------------------------------------------------------\r\n');
for n = 1:37
    line = fgetl(fid1);
    fprintf(fid2, '%s\r\n', line);
end
fclose(fid1);
fclose(fid2);

end % function writeInpFileNewMain