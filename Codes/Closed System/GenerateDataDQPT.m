function DATA = GenerateDataDQPT(J,h,alpha,Nspins)


%% Nalpha: factor used for the renormalization
Nalpha = 0;
for i=1:Nspins
    for j=1:Nspins
        if i~=j
            Nalpha = Nalpha + abs(i-j)^(-alpha)/(Nspins-1);
        end
    end
end

%% Spin operators for S=1/2
Sx = [0 1;1 0];
Sy = [0 -1i;1i 0];
Sz = [1 0;0 -1];

%% Hamiltonian H0
H0 = 0;
for i=1:Nspins
    Sxi = getSci(Sx,i,Nspins,2);
    for j=1:Nspins
        if i~=j
            Sxj = getSci(Sx,j,Nspins,2);
            H0 = H0 - J*abs(i-j)^(-alpha)/Nalpha*Sxi*Sxj; % Interaction Hamiltonian 
        end
    end
end

%% Spins matrices for each site
SSx = cell(1,Nspins);
SSy = cell(1,Nspins);
SSz = cell(1,Nspins);
for i=1:Nspins
    SSx{i} = getSci(Sx,i,Nspins,2);
    SSy{i} = getSci(Sy,i,Nspins,2);
    SSz{i} = getSci(Sz,i,Nspins,2);
end

%% Eigenstates Sx
xr = [1 1]'/sqrt(2);        % Single-particle state |Psi_{-->}>
xl = [-1 1]'/sqrt(2);       % Single-particle state |Psi_{<--}>

%% Eigenstates and eigenvalues of H0 (spin-spin interaction)
Xr = xr;
Xl = xl;
for n=1:Nspins-1
    Xr = kron(Xr,xr);     % Many-body state |Psi_{-->}>
    Xl = kron(Xl,xl);     % Many-body state |Psi_{<--}>
end

%% Initial condition
PSI_0 = Xl;       

%% Spherical angles used as data
N = 200;
phi = linspace(0,pi,N);
theta = linspace(0,2*pi,N);

%% Tuplas used to make faster the algorithm
[ca, cb] = ndgrid(theta,phi);
combs = [ca(:), cb(:)]; 
theta = combs(:,1); phi = combs(:,2);
DATA = zeros(length(combs),3);
parfor i=1:length(combs)
    DATA(i,:)= [theta(i) phi(i) getDQPT_SystemDynamics(h,theta(i),phi(i),J,Nspins,SSx, SSy, SSz, PSI_0, Xr, Xl, H0)];
end

%% Figure to obser the pattern of dynamical singularities
figure()
box on
gscatter(DATA(:,1)/pi,DATA(:,2)/pi,DATA(:,3));
title('DQPT phase diagram')
xlabel('$\theta/\pi$','Interpreter','latex','Fontsize', 21)
ylabel('$\phi/\pi$','Interpreter','latex','Fontsize', 21)
set(gca,'fontsize',21)

end