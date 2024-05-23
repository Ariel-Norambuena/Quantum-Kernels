function DQPT = getDQPT_MarkovianNoise(B,theta,phi,alpha,J,Nspins,Nalpha,Sx,Sy,Sz,xr,xl)

%% Magnetic field
Bx = B*cos(theta)*sin(phi);
By = B*sin(theta)*sin(phi);
Bz = B*cos(phi);

%% Hamiltonians H0 and H1
H0 = 0; H1 = 0;
for i=1:Nspins
    Sxi = getSci(Sx,i,Nspins,2);
    Syi = getSci(Sy,i,Nspins,2);
    Szi = getSci(Sz,i,Nspins,2);
    H1 = H1-Bx*Sxi-By*Syi-Bz*Szi;       % Hamiltonian for the magnetic field
    for j=1:Nspins
        if i~=j
            Sxj = getSci(Sx,j,Nspins,2);
            H0 = H0 - J*abs(i-j)^(-alpha)/Nalpha*Sxi*Sxj; % Interaction Hamiltonian
        end
    end
end

%% Eigenstates and eigenvalues of H0 (spin-spin interaction)
Xr = xr;
Xl = xl;
for n=1:Nspins-1
    Xr = kron(Xr,xr);        % Many-body state |Psi_{-->}>
    Xl = kron(Xl,xl);        % Many-body state |Psi_{<--}>
end

%% First: we consider the ground state wavefunction
PSI_0 = Xl;         % Initial condition

%% Time evolution
Nt = 1000;
ti = 1e-2;                  % Initial time
if B ==0
    wmin = J;
else
    wmin = min(J,B);
end
tf = 20*pi/wmin;
dt = (tf-ti)/(Nt-1);        % Step time
t = ti:dt:tf;               % Time vector

%% Populations
Pr = zeros(size(t));
Pl = zeros(size(t));

%% Hamiltonian
H = H0 + H1 ;                           % Total hamiltonian
Is = eye(2^Nspins);
L_H = -1i*kron(Is,H)+1i*kron(H.',Is);

%% Losses for amplitude damping
gamma = 0.02*J;     % decay rate
dimT = 2^Nspins;    % total dimension Hilbert space 
L_diss = 0;
for i=1:Nspins
    L_diss = 0;
    Sxi = getSci(Sx,i,Nspins,2);
    Syi = getSci(Sy,i,Nspins,2);
    Smi = Sxi-1i*Syi;
    L_diss = L_diss + LindbladSuperOperator(gamma,Smi,dimT);
end
L = L_H + L_diss;
sigma_rr = Xr*Xr';
sigma_ll = Xl*Xl';
rho_0 = PSI_0*PSI_0';
rhot = SolvesMasterEquation(L,rho_0,t);

for n=1:length(t)
    rho = rhot{n};
    Pr(n) = real(trace(rho*sigma_rr));
    Pl(n) = real(trace(rho*sigma_ll));
end

%% Determine if exist DQPT or not
Diff = Pl-Pr;
if ~isempty(find(Diff<0, 1))
    DQPT= 1;  % Yes
else
    DQPT= -1; % Not
end

end