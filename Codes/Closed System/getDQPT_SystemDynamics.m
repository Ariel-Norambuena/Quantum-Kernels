function DQPT = getDQPT_SystemDynamics(B,theta,phi,J,Nspins, SSx, SSy, SSz, PSI_0, Xr, Xl,H0)

%% Magnetic field 
Bx = B*cos(theta)*sin(phi);
By = B*sin(theta)*sin(phi);
Bz = B*cos(phi);

%% Hamiltonians H0 and H1
H1 = 0;
for i=1:Nspins
    Sxi = SSx{i};
    Syi = SSy{i};
    Szi = SSz{i};
    H1 = H1-Bx*Sxi-By*Syi-Bz*Szi;       % Hamiltonian for the magnetic field
end

%% Time evolution
H = H0 + H1;                % Total hamiltonian
Nt = 1000;                  % Number steps for time
ti = 1e-2;                  % Initial time
if B ==0
    wmin = J;
else
    wmin = min(J,B);
end
tf = 20*pi/wmin;            % Large time window to explore DQPT
dt = (tf-ti)/(Nt-1);        % Step time
t = ti:dt:tf;               % Time vector
U = expm(-1i*H*dt);         % Time propagator operator U(dt)

%% Variables to save information of the return probabilities
Pr = zeros(size(t));
Pl = zeros(size(t));

for n=1:length(t)
    if n==1
        PSI = PSI_0;        % Initial wave function
    else
        PSI = U*PSI;        % Wave function at time t_n = n*dt + t0
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