function writeInpFileBModes(SIM, ...
                            ANLS, ...
                            BLADE, ...
                            elas_tw, ...
                            iner_tw, ...
                            mass_den, ...
                            flapIner_cm, ...
                            edgeIner_cm, ...
                            flapEI_sc, ...
                            edgeEI_sc, ...
                            tor_stff, ...
                            axial_stff, ...
                            cm_offst, ...
                            sc_offst, ...
                            tc_offst)

%% re-assign some structure variable names (for convenience)                       
N_MODES    = ANLS.N_MODES;
N_ELEMS    = ANLS.N_ELEMS;
NUM_SEC    = BLADE.NUM_SEC;
zFrac      = BLADE.zFrac;
ROT_SPD    = BLADE.ROT_SPD;
HUB_RAD    = BLADE.HUB_RAD;
BLD_LENGTH = BLADE.BLD_LENGTH;
PRE_CONE   = BLADE.PRE_CONE;
BLD_PITCH  = BLADE.BLD_PITCH;

%%                  
date = datestr(now, 'mmmm dd, yyyy HH:MM AM');
                        
%% write the BModes section properties file
fid = fopen([SIM.outputDir filesep SIM.case '_BModesProps.bmi'], 'w');                        
   
fprintf(fid,'Blade section properties.  Generated on %s by Co-Blade v%s\r\n', date, SIM.version);
fprintf(fid,'%-9.0f  n_secs:     number of blade sections at which properties are specified (-)\r\n', NUM_SEC);
fprintf(fid,'\r\n');
fprintf(fid,'sec_loc  str_tw  tw_iner  mass_den flp_iner  edge_iner  flp_stff    edge_stff     tor_stff    axial_stff  cg_offst  sc_offst tc_offst\r\n');
fprintf(fid,'(-)      (deg)    (deg)   (kg/m)    (kg-m)    (kg-m)     (Nm^2)      (Nm^2)        (Nm^2)        (N)         (m)       (m)     (m)\r\n');
for i = 1:NUM_SEC
    fprintf(fid,'%8.3f  %8.3f  %8.3f  %8.3E  %8.3E  %8.3E  %8.3E  %8.3E  %8.3E  %8.3E  %8.4f  %8.4f  %8.4f\r\n', ...  
                 zFrac(i),elas_tw(i),iner_tw(i),mass_den(i),flapIner_cm(i),edgeIner_cm(i),flapEI_sc(i),edgeEI_sc(i),tor_stff(i),axial_stff(i),cm_offst(i),sc_offst(i),tc_offst(i));
end

fclose(fid);
             
%% write the BModes main input file
fid = fopen([SIM.outputDir filesep SIM.case '_BModes.bmi'], 'w');  

fprintf(fid,'======================   BModes v3.00 Main Input File  ==================\r\n');
fprintf(fid,'Generated on %s by Co-Blade v%s \r\n', date, SIM.version);
fprintf(fid,'\r\n');
fprintf(fid,'--------- General parameters ---------------------------------------------------------------------\r\n');
fprintf(fid,'false     Echo        Echo input file contents to *.echo file if true.\r\n');
fprintf(fid,'1         beam_type   1: blade, 2: tower (-)\r\n');
fprintf(fid,'%-9.3f rot_rpm:    rotor speed, automatically set to zero for tower modal analysis (rpm)\r\n', ROT_SPD);
fprintf(fid,'1.0       rpm_mult:   rotor speed muliplicative factor (-)\r\n');
fprintf(fid,'%-9.3f radius:     rotor tip radius measured along coned blade axis OR tower height (m)\r\n', HUB_RAD + BLD_LENGTH);
fprintf(fid,'%-9.3f hub_rad:    hub radius measured along coned blade axis OR tower rigid-base height (m)\r\n', HUB_RAD);
fprintf(fid,'%-9.3f precone:    built-in precone angle, automatically set to zero for a tower (deg)\r\n', PRE_CONE);
fprintf(fid,'%-9.3f bl_thp:     blade pitch setting, automatically set to zero for a tower (deg)\r\n', BLD_PITCH);
fprintf(fid,'1         hub_conn:   hub-to-blade or tower-base boundary condition [1: cantilevered; 2: free-free; 3: only axial and torsion constraints] (-)\r\n');
fprintf(fid,'%-9.0f modepr:     number of modes to be printed (-)\r\n', N_MODES);
fprintf(fid,'true      TabDelim    (true: tab-delimited output tables; false: space-delimited tables)\r\n');
fprintf(fid,'false     mid_node_tw (true: output twist at mid-node of elements; false: no mid-node outputs)\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'--------- Blade-tip or tower-top mass properties --------------------------------------------\r\n');
fprintf(fid,'0.0       tip_mass    blade-tip or tower-top mass (kg)\r\n');
fprintf(fid,'0.0       cm_loc      tip-mass c.m. offset from the blade axis measured along the tip section y reference axis (m)\r\n');
fprintf(fid,'0.0       cm_axial    tip-mass c.m. offset tower tip measures axially along the z axis (m)\r\n');
fprintf(fid,'0.0       ixx_tip     blade lag mass moment of inertia about the tip-section x reference axis (kg-m^2)\r\n');
fprintf(fid,'0.0       iyy_tip     blade flap mass moment of inertia about the tip-section y reference axis (kg-m^2)\r\n');
fprintf(fid,'0.0       izz_tip     torsion mass moment of inertia about the tip-section z reference axis (kg-m^2)\r\n');
fprintf(fid,'0.0       ixy_tip     cross product of inertia about x and y reference axes(kg-m^2)\r\n');
fprintf(fid,'0.0       izx_tip     cross product of inertia about z and x reference axes(kg-m^2)\r\n');
fprintf(fid,'0.0       iyz_tip     cross product of inertia about y and z reference axes(kg-m^2)\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'--------- Distributed-property identifiers --------------------------------------------------------\r\n');
fprintf(fid,'1         id_mat:     material_type [1: isotropic; non-isotropic composites option not yet available]\r\n');
fprintf(fid,'''%s'' sec_props_file   name of beam section properties file (-)\r\n', [SIM.case '_BModesProps.bmi']);
fprintf(fid,'\r\n');
fprintf(fid,'Property scaling factors..............................\r\n');
fprintf(fid,'1.0       sec_mass_mult:   mass density multiplier (-)\r\n');
fprintf(fid,'1.0       flp_iner_mult:   blade flap or tower f-a inertia multiplier (-)\r\n');
fprintf(fid,'1.0       lag_iner_mult:   blade lag or tower s-s inertia multiplier (-)\r\n');
fprintf(fid,'1.0       flp_stff_mult:   blade flap or tower f-a bending stiffness multiplier (-)\r\n');
fprintf(fid,'1.0       edge_stff_mult:  blade lag or tower s-s bending stiffness multiplier (-)\r\n');
fprintf(fid,'1.0       tor_stff_mult:   torsion stiffness multiplier (-)\r\n');
fprintf(fid,'1.0       axial_stff_mult: axial stiffness multiplier (-)\r\n');
fprintf(fid,'1.0       cg_offst_mult:   cg offset multiplier (-)\r\n');
fprintf(fid,'1.0       sc_offst_mult:   shear center multiplier (-)\r\n');
fprintf(fid,'1.0       tc_offst_mult:   tension center multiplier (-)\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'--------- Finite element discretization --------------------------------------------------\r\n');
fprintf(fid,'%-9.0f nselt:     no of blade or tower elements (-)\r\n', N_ELEMS);
fprintf(fid,'Distance of element boundary nodes from blade or flexible-tower root (normalized wrt blade or tower length), el_loc()\r\n');
fprintf(fid,'%-9.3f', linspace(0, 1, N_ELEMS + 1));
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'END of Main Input File Data *********************************************************************\r\n');
fprintf(fid,'*************************************************************************************************\r\n');

fclose(fid);

end % function writeInpFileBModes

