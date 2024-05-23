function DQPT = getDQPT_OscillatingControl(B,theta,phi,alpha,J,Nspins,Nalpha,Sx,Sy,Sz,xr,xl,SSx,SSy,SSz)

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
ti = 1e-2;                     % Initial time
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

%% Time-dependent magnetic field
w = 1.2*J;              % frequency
Bz0 = 0.05*B;           % amplitude
B_osc_z = Bz0*sin(w*t); % time-dependent profile

for n=1:length(t)
    if n==1
        PSI = PSI_0;        % Initial wavefunction
    else
        Hosc  = B_osc_z(n)*SSz;                           % Hamiltonian for the magnetic field
        H = H0 + H1 + Hosc;                           % Total hamiltonian
        U = expm(-1i*H*dt);                             % Time propagator operator U(dt)
        PSI = U*PSI;        % Wavefunction at time t_n
    end
    Pr(n) = real(abs(Xr'*PSI)^2);
    Pl(n) = real(abs(Xl'*PSI)^2);
end

%% Determine if exist DQPT or not
Diff = Pl-Pr;
if ~isempty(find(Diff<0, 1))
    DQPT= 1;  % Yes
else
    DQPT= -1; % Not
end

end