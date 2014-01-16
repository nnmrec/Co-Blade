function CLT = cltAnalysis(lamNum,nPlies,tPly,fibAng,E11,E22,G12,nu12,density)
    
d2r = pi/180;   % convert degrees to radians

%% Compute total thickness of the laminate
t = sum( nPlies.*tPly );

%% Compute mass per unit volume
massByV = sum(nPlies.*tPly.*density) / t; % (kg/m^3)

%% Compute ply (lamina) interface locations
N    = length(lamNum);
z    = zeros(N+1,1);
z(1) = -t/2;
for k = 1:N
    z(k+1) = z(k) + nPlies(k)*tPly(k);
end

%% Calculate the reduced stiffness matrix, Q, for each lamina
Q = zeros(3,3,N);
for k = 1:N    
    Q(1,1,k) = E11(k)^2              / (E11(k) - nu12(k)^2 * E22(k));
    Q(2,2,k) = E11(k)*E22(k)         / (E11(k) - nu12(k)^2 * E22(k));
    Q(1,2,k) = nu12(k)*E11(k)*E22(k) / (E11(k) - nu12(k)^2 * E22(k));
    Q(2,1,k) = Q(1,2,k);
    Q(3,3,k) = G12(k);     
end 

%% Calculate the transformed reduced stiffness matrix, Qbar, and transformation matrix, T, for each lamina
Qbar = zeros(3,3,N);
T    = zeros(3,3,N);
for k = 1:N
    c = cos( fibAng(k)*d2r );
    s = sin( fibAng(k)*d2r );
    
    Qbar(1,1,k) = (c^4)*Q(1,1,k) + 2*(c^2)*(s^2)*(Q(1,2,k) + 2*Q(3,3,k)) + (s^4)*Q(2,2,k);
    Qbar(1,2,k) = (c^4 + s^4)*Q(1,2,k) + (c^2)*(s^2)*(Q(1,1,k) + Q(2,2,k) - 4*Q(3,3,k));
    Qbar(2,1,k) = Qbar(1,2,k);
    Qbar(1,3,k) = (c^3)*s*(Q(1,1,k) - Q(1,2,k) - 2*Q(3,3,k)) - c*(s^3)*(Q(2,2,k) - Q(1,2,k) - 2*Q(3,3,k));
    Qbar(3,1,k) = Qbar(1,3,k);
    Qbar(2,2,k) = (s^4)*Q(1,1,k) + 2*(c^2)*(s^2)*(Q(1,2,k) + 2*Q(3,3,k)) + (c^4)*Q(2,2,k);
    Qbar(2,3,k) = c*(s^3)*(Q(1,1,k) - Q(1,2,k) - 2*Q(3,3,k)) - s*(c^3)*(Q(2,2,k) - Q(1,2,k) - 2*Q(3,3,k));
    Qbar(3,2,k) = Qbar(2,3,k);
    Qbar(3,3,k) = (c^2)*(s^2)*(Q(1,1,k) + Q(2,2,k) - 2*Q(1,2,k) - 2*Q(3,3,k)) + (c^4 + s^4)*Q(3,3,k);
    
    T(1,1,k) = c^2;
    T(1,2,k) = s^2;
    T(1,3,k) = 2*c*s;
    T(2,1,k) = s^2;
    T(2,2,k) = c^2;
    T(2,3,k) = -2*c*s;
    T(3,1,k) = -c*s;
    T(3,2,k) = c*s;
    T(3,3,k) = c^2 - s^2;
 
end 

%% Calculate the A, B, D, and ABD matrices
A = zeros(3,3); % extensional stiffness matrix
B = zeros(3,3); % coupling matrix
D = zeros(3,3); % bending stiffness matrix
for i = 1:3
    for j = 1:3
        sumA = 0;
        sumB = 0;
        sumD = 0;
        for k = 1:N
            sumA = sumA + Qbar(i,j,k)*(z(k+1) - z(k));
            sumB = sumB + Qbar(i,j,k)*(z(k+1)^2 - z(k)^2);
            sumD = sumD + Qbar(i,j,k)*(z(k+1)^3 - z(k)^3);
        end
        A(i,j) = sumA;
        B(i,j) = (1/2) * sumB;
        D(i,j) = (1/3) * sumD;     
    end
end
ABD = [A B; 
       B D];

%% The effective elastic properties of the laminate
a      = inv(A);
E_eff  = 1/(t*a(1,1));   % effective Young's modulus of the laminate
G_eff  = 1/(t*a(3,3));   % effective shear modulus of the laminate
nu_eff = -a(1,2)/a(1,1); % effective Poisson ratio of the laminate

%% Collect output
CLT.t       = t;
CLT.E_eff   = E_eff;
CLT.G_eff   = G_eff;
CLT.nu_eff  = nu_eff;
CLT.massByV = massByV;
CLT.A       = A;
CLT.B       = B;
CLT.D       = D;
CLT.ABD     = ABD;
CLT.Qbar    = Qbar;
CLT.T       = T;
CLT.z_i     = z;

end %function cltAnalysis

