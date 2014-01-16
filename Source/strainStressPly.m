function Out = strainStressPly(Panel, MidPlane, MATS)

%% re-assign some structure variable names (for convenience)
s11_fT  = MATS.s11_fT;
s11_fC  = MATS.s11_fC;
s22_yT  = MATS.s22_yT;
s22_yC  = MATS.s22_yC;
s12_y   = MATS.s12_y;
nPanels = Panel.nPanels;

%%
if nPanels < 1
    Out.e_11_maxabs = [];
    Out.e_22_maxabs = [];
    Out.e_12_maxabs = [];
    Out.s_11_maxabs = [];
    Out.s_22_maxabs = [];
    Out.s_12_maxabs = [];
    Out.s_11_fc_T   = [];
    Out.s_11_fc_C   = [];
    Out.s_22_fc_T   = [];
    Out.s_22_fc_C   = [];
    Out.s_12_fc_S   = [];
    return
end

Out(nPanels, 1) = struct('e_11_maxabs', [], ...
                         'e_22_maxabs', [], ...
                         'e_12_maxabs', [], ...
                         's_11_maxabs', [], ...
                         's_22_maxabs', [], ...
                         's_12_maxabs', [], ...
                         's_11_fc_T', [], ...
                         's_11_fc_C', [], ...
                         's_22_fc_T', [], ...
                         's_22_fc_C', [], ...
                         's_12_fc_S', []);
for n = 1:nPanels

    z_i  = Panel.z_i{n};     % plate z-coordinate of lamina interface locations
    nZ   = length(z_i);      % number of lamina interfaces in the panel (laminate)
    nLam = nZ - 1;           % number of laminas in the panel (laminate)
    
    % calculate the strains at the lamina interface locations in the panel (laminate)
    e_xx = cell(nZ, 1); % normal strain in the x-direction at the i-th lamina interface
    e_yy = cell(nZ, 1); % normal strain in the y-direction at the i-th lamina interface
    e_xy = cell(nZ, 1); % shear strain w.r.t. x-y at the i-th lamina interface
    for i = 1:nZ
        e_xx{i} = MidPlane.eo_xx{n} + z_i(i)*MidPlane.k_xx{n};  
        e_yy{i} = MidPlane.eo_yy{n} + z_i(i)*MidPlane.k_yy{n};  
        e_xy{i} = MidPlane.eo_xy{n} + z_i(i)*MidPlane.k_xy{n};  
    end
    
    e_11_maxabs = cell(nLam, 1); % 1st principal  strain in the i-th lamina which has the largest absolute value
    e_22_maxabs = cell(nLam, 1); % 2nd principal strain in the i-th lamina which has the largest absolute value
    e_12_maxabs = cell(nLam, 1); % maximum absolute value of principal shear strain in the i-th lamina
    s_11_maxabs = cell(nLam, 1); % 1st principal  stress in the i-th lamina which has the largest absolute value
    s_22_maxabs = cell(nLam, 1); % 2nd principal stress in the i-th lamina which has the largest absolute value
    s_12_maxabs = cell(nLam, 1); % maximum absolute value of principal shear stress in the i-th lamina which
    s_11_fc_T   = cell(nLam, 1); % failure criteria for 1st principal stress (failure in tension)
    s_11_fc_C   = cell(nLam, 1); % failure criteria for 1st principal stress (failure in compression)
    s_22_fc_T   = cell(nLam, 1); % failure criteria for 2nd principal stress (yield in tension)
    s_22_fc_C   = cell(nLam, 1); % failure criteria for 2nd principal stress (yield in compression)
    s_12_fc_S   = cell(nLam, 1); % failure criteria for principal shear stress (yield in shear)
    for i = 1:nLam
        % in each lamina, convert the strains into principal strains
        e_p_ou  = Panel.T{n}(:,:,i) * [e_xx{i}';   e_yy{i}';   e_xy{i}'./2];    % principal strains on the outer surface of the i-th lamina
        e_p_in  = Panel.T{n}(:,:,i) * [e_xx{i+1}'; e_yy{i+1}'; e_xy{i+1}'./2];  % principal strains on the inner surface of the i-th lamina
        e_11_ou = e_p_ou(1,:)';     %  first principal strain on the outer surface of the i-th lamina
        e_22_ou = e_p_ou(2,:)';     % second principal strain on the outer surface of the i-th lamina
        e_12_ou = 2.*e_p_ou(3,:)';  %  principal shear strain on the outer surface of the i-th lamina
        e_11_in = e_p_in(1,:)';     %  first principal strain on the inner surface of the i-th lamina
        e_22_in = e_p_in(2,:)';     % second principal strain on the inner surface of the i-th lamina
        e_12_in = 2.*e_p_in(3,:)';  %  principal shear strain on the inner surface of the i-th lamina
        % principal strains
        e_11_maxabs{i} = maxAbs(e_11_ou, e_11_in);
        e_22_maxabs{i} = maxAbs(e_22_ou, e_22_in);
        e_12_maxabs{i} = max(abs(e_12_ou), abs(e_12_in)); 
        
        % calculate the stress in each lamina, at the lamina interface locations
        s_ou    = Panel.Qbar{n}(:,:,i) * [e_xx{i}';   e_yy{i}';   e_xy{i}'];   % stresses on the outer surface of the i-th lamina
        s_in    = Panel.Qbar{n}(:,:,i) * [e_xx{i+1}'; e_yy{i+1}'; e_xy{i+1}']; % stresses on the inner surface of the i-th lamina
        s_xx_ou = s_ou(1,:)'; % normal stresses in the x-direction on the outer surface of the i-th lamina
        s_yy_ou = s_ou(2,:)'; % normal stresses in the y-direction on the outer surface of the i-th lamina
        s_xy_ou = s_ou(3,:)'; % shear stresses w.r.t. x-y on the outer surface of the i-th lamina
        s_xx_in = s_in(1,:)'; % normal stresses in the x-direction on the inner surface of the i-th lamina
        s_yy_in = s_in(2,:)'; % normal stresses in the y-direction on the inner surface of the i-th lamina
        s_xy_in = s_in(3,:)'; % shear stresses w.r.t. x-y on the inner surface of the i-th lamina
            
        % in each lamina convert the stresses into principle stresses
        s_p_ou  = Panel.T{n}(:,:,i) * [s_xx_ou';   s_yy_ou';   s_xy_ou'];  % principal stresses on the outer surface of the i-th lamina
        s_p_in  = Panel.T{n}(:,:,i) * [s_xx_in';   s_yy_in';   s_xy_in'];  % principal stresses on the inner surface of the i-th lamina
        s_11_ou = s_p_ou(1,:)';   % first principal stress  on the outer surface of the i-th lamina
        s_22_ou = s_p_ou(2,:)';   % second principal stress on the outer surface of the i-th lamina
        s_12_ou = s_p_ou(3,:)';   % principal shear stress  on the outer surface of the i-th lamina
        s_11_in = s_p_in(1,:)';   % first principal stress  on the inner surface of the i-th lamina
        s_22_in = s_p_in(2,:)';   % second principal stress on the inner surface of the i-th lamina
        s_12_in = s_p_in(3,:)';   % principal shear stress  on the inner surface of the i-th lamina
        % principal stresses
        s_11_maxabs{i} = maxAbs(s_11_ou, s_11_in);
        s_22_maxabs{i} = maxAbs(s_22_ou, s_22_in);
        s_12_maxabs{i} = max(abs(s_12_ou), abs(s_12_in));
        
        % stress failure criteria
        s_11_fc_T{i} = max(0, s_11_maxabs{i}) ./ s11_fT(Panel.matID{n}(i));
        s_11_fc_C{i} = min(0, s_11_maxabs{i}) ./ s11_fC(Panel.matID{n}(i));
        s_22_fc_T{i} = max(0, s_22_maxabs{i}) ./ s22_yT(Panel.matID{n}(i));
        s_22_fc_C{i} = min(0, s_22_maxabs{i}) ./ s22_yC(Panel.matID{n}(i));
        s_12_fc_S{i} = s_12_maxabs{i} ./ s12_y(Panel.matID{n}(i));
        
    end
    
    %% Collect the output
    Out(n).e_11_maxabs = e_11_maxabs;
    Out(n).e_22_maxabs = e_22_maxabs;
    Out(n).e_12_maxabs = e_12_maxabs;    
    Out(n).s_11_maxabs = s_11_maxabs;
    Out(n).s_22_maxabs = s_22_maxabs;
    Out(n).s_12_maxabs = s_12_maxabs;
    Out(n).s_11_fc_T   = s_11_fc_T;
    Out(n).s_11_fc_C   = s_11_fc_C;
    Out(n).s_22_fc_T   = s_22_fc_T;
    Out(n).s_22_fc_C   = s_22_fc_C;
    Out(n).s_12_fc_S   = s_12_fc_S;


end 

end % function strainStressPly

function f = maxAbs(a, b)

N = numel(a);
f = zeros(N, 1);
for i = 1:N   
   if abs(a(i)) > abs(b(i))
       f(i) = a(i);
   else
       f(i) = b(i);
   end
end

end % function maxAbs