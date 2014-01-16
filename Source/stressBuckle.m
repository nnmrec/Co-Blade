function Buckle = stressBuckle(Panel, NormS, ShearS, flag)

nPanels = Panel.nPanels;

if nPanels < 1
    Buckle = [];
    return
end

switch flag
    case {'top','bottom'}
        
        avg_stress_zz = zeros(nPanels, 1);
        avg_stress_zs = zeros(nPanels, 1);
        avg_r         = zeros(nPanels, 1);
        b             = zeros(nPanels, 1);
        for n = 1:nPanels
            avg_stress_zz(n) = mean(NormS.stress_zz{n});
            avg_stress_zs(n) = mean(abs(ShearS.stress_zs{n}));
            avg_r(n)         = radiusCurvature(Panel.xs{n}, Panel.ys{n});
            b(n)             = Panel.s{n}(end);
        end
        
        % overwrite any positive values, we only care about compressive stresses here
        avg_stress_zz(avg_stress_zz > 0) = 0;
        
        t      = Panel.t;
        E_eff  = Panel.E_eff;
        nu_eff = Panel.nu_eff;
        sroot  = sqrt( 12.*(1-nu_eff.^2).*((t./avg_r).^2) + (pi.*t./b).^4 );
        
        % critical buckling compressive normal stress
        stress_zz_crit = -(1/6) .* E_eff./(1-nu_eff.^2) .* (sroot + (pi.*t./b).^2);
        
        % critical buckling shear stress
        stress_zs_crit = 0.1 .* E_eff .* (t./avg_r) + 5 .* E_eff .* ((t./b).^2);
        
        R_comp  = avg_stress_zz ./ stress_zz_crit;
        R_shear = avg_stress_zs ./ stress_zs_crit;
        Buckle  = R_comp + R_shear.^(3/2);
        
    case 'web'
        
        min_stress_zz = zeros(nPanels, 1); % maximum compressive stress from bending
        avg_stress_zs = zeros(nPanels, 1);
        b             = zeros(nPanels, 1);
        for n = 1:nPanels
            min_stress_zz(n) = min(NormS.stress_zz{n});
            avg_stress_zs(n) = mean(abs(ShearS.stress_zs{n}));
            b(n)             = Panel.s{n}(end);
        end
        
        % overwrite any positive values, we only care about compressive stresses here
        min_stress_zz(min_stress_zz > 0) = 0;
        
        t      = Panel.t;
        E_eff  = Panel.E_eff;
        nu_eff = Panel.nu_eff;
        
        % critical buckling compressive normal stress
        stress_zz_crit = -21.75 .* E_eff .* (t./b).^2;
        
        % critical buckling shear stress
        stress_zs_crit = 4.4 .* E_eff./(1-nu_eff.^2) .* (t./b).^2;
        
        R_bend  = min_stress_zz ./ stress_zz_crit;
        R_shear = avg_stress_zs ./ stress_zs_crit;
        Buckle  = R_bend.^2 + R_shear.^2;
        
end

end % function stressBuckle
