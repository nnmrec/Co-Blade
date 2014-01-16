function Out = strainPlate(Panel, NormS, ShearS)

nPanels = Panel.nPanels;

if nPanels < 1                        
    Out.eo_xx = [];
    Out.eo_yy = []; 
    Out.eo_xy = [];
    Out.k_xx  = [];
    Out.k_yy  = []; 
    Out.k_xy  = [];   
    return
end

eo_xx = cell(nPanels, 1);  % plate midplane normal strain in the x-direction
eo_yy = cell(nPanels, 1);  % plate midplane normal strain in the y-direction
eo_xy = cell(nPanels, 1);  % plate midplane shear strain w.r.t. the x-y axes
k_xx  = cell(nPanels, 1);  % plate midplane curvature w.r.t. the x-axis
k_yy  = cell(nPanels, 1);  % plate midplane curvature w.r.t. the y-axis
k_xy  = cell(nPanels, 1);  % plate midplane curvature w.r.t. the x-y axes
for n = 1:nPanels
    t    = Panel.t(n);
    ABD  = Panel.ABD{n};
    s_zz = NormS.stress_zz{n};  % normal stresses computed by beam theory around the perimeter of the panel
    s_zs = ShearS.stress_zs{n}; % shear stress computed by beam theory at the midwall of the panel
    
    nS      = length(s_zz) / 2;          % the number of points along the midwall where midplane strains and curvatures will be calculated
    s_zz_ou = s_zz(1:nS);                % normal stress on the outside blade surface from beam theory (N/m^2)
    s_zz_in = flipud( s_zz(nS+1:end) );  % normal stress on the inside blade surface from beam theory (N/m^2)
    
    % equivalent distributed forces and moments applied along the plate edge
    Nxx = t.*(s_zz_ou + s_zz_in)./2;                % in-plane distributed normal force per unit length (N/m) in the x-direction
    Nyy = 0;                                        % in-plane distributed normal force per unit length (N/m) in the y-direction
    Nxy = t.*s_zs;                                  % in-plane distributed shear force per unit length (N/m) w.r.t. x-y axes
    Mxx = t.*(s_zz_ou.*(-t./2) + s_zz_in.*t./2)./2; % distributed bending moment per unit length (Nm/m) w.r.t. the plate y axis 
    Myy = 0;                                        % distributed bending moment per unit length (Nm/m) w.r.t. the plate x axis 
    Mxy = 0;                                        % distributed in-plane torque moment per unit length (Nm/m)
    
    % plate midplane strains and curvatures
    eo_xx{n} = zeros(nS, 1);  % plate midplane normal strain in the x-direction
    eo_yy{n} = zeros(nS, 1);  % plate midplane normal strain in the y-direction
    eo_xy{n} = zeros(nS, 1);  % plate midplane shear strain w.r.t. the x-y axes
    k_xx{n}  = zeros(nS, 1);  % plate midplane curvature w.r.t. the x-axis
    k_yy{n}  = zeros(nS, 1);  % plate midplane curvature w.r.t. the y-axis
    k_xy{n}  = zeros(nS, 1);  % plate midplane curvature w.r.t. the x-y axes
    for i = 1:nS
        NM       = [Nxx(i); Nyy; Nxy(i); Mxx(i); Myy; Mxy];
        ek       = ABD \ NM;
        eo_xx{n}(i) = ek(1);
        eo_yy{n}(i) = ek(2);
        eo_xy{n}(i) = ek(3);
        k_xx{n}(i)  = ek(4);
        k_yy{n}(i)  = ek(5);
        k_xy{n}(i)  = ek(6);
    end
end

%% Collect the output
Out.eo_xx = eo_xx;
Out.eo_yy = eo_yy; 
Out.eo_xy = eo_xy;
Out.k_xx  = k_xx;
Out.k_yy  = k_yy; 
Out.k_xy  = k_xy;
    
end % function strainPlate

