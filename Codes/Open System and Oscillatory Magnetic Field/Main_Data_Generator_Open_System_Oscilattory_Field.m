tic
clear all 
close all

%% Physical parameters
J = 1;
alpha = 0.5;
h = 0.95*J;
Nspins = 2;

%% N_alpha: renormalization factor
Nalpha = 0;
for i=1:Nspins
    for j=1:Nspins
        if i~=j
            Nalpha = Nalpha + abs(i-j)^(-alpha)/(Nspins-1);
        end
    end
end

%% Spin operators
Sx = [0 1;1 0];
Sy = [0 -1i;1i 0];
Sz = [1 0;0 -1];

%% Noise operator
SSz = 0;SSy = 0;SSx = 0;
for i=1:Nspins
    Sxi = getSci(Sx,i,Nspins,2);
    Syi = getSci(Sy,i,Nspins,2);
    Szi = getSci(Sz,i,Nspins,2);
    SSx  = SSx + Sxi;       % Hamiltonian for the magnetic field
    SSy  = SSy + Syi;       % Hamiltonian for the magnetic field
    SSz  = SSz + Szi;       % Hamiltonian for the magnetic field
end

%% Eigenstates Sx
xr = [1 1]'/sqrt(2);        % Single-particle state |Psi_{-->}>
xl = [-1 1]'/sqrt(2);       % Single-particle state |Psi_{<--}>

%% Spherical angles as data
N = 100;
theta = linspace(0,2*pi,N);
phi = linspace(0,pi,N);

B0 = linspace(0,J,N);
w = linspace(0,J,N);

% Tuplas for O(N) complexity: faster than two for
[A,B] = meshgrid(theta,phi);
c=cat(2,A',B');
d=reshape(c,[],2);
theta = d(:,1); phi = d(:,2);

[A,B] = meshgrid(B0,w);
c=cat(2,A',B');
d=reshape(c,[],2);
B0 = d(:,1); w = d(:,2);

DATA1 = zeros(N*N,3);
DATA2 = zeros(N*N,3);

parfor i=1:N*N
    DATA1(i,:) = [theta(i) phi(i) getDQPT_OscillatingControl(h,theta(i),phi(i),alpha,J,Nspins,Nalpha,Sx,Sy,Sz,xr,xl,SSx,SSy,SSz)];
    DATA2(i,:) = [theta(i) phi(i) getDQPT_MarkovianNoise(h,theta(i),phi(i),alpha,J,Nspins,Nalpha,Sx,Sy,Sz,xr,xl)];
end

writematrix(DATA1,'DATA_DQPT_N=2_theta_phi_OscillatoryField.txt')
writematrix(DATA2,'DATA_DQPT_N=2_theta_phi_MarkovianNoise.txt')

opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 8;
opts.height     = 6;
opts.fontType   = 'Times';
opts.fontSize   = 9;

figure()
box on
s1 = gscatter(DATA1(:,1)/pi,DATA1(:,2)/pi,DATA1(:,3));
s1(1).MarkerSize = 7; s1(1).Color = 'red';
s1(2).MarkerSize = 7; s1(2).Color = 'blue';
title('$\mbox{Oscillating field}$','Interpreter','latex','Fontsize', 15)
xlabel('$\theta_m/\pi$','Interpreter','latex','Fontsize', 15)
ylabel('$\phi_m/\pi$','Interpreter','latex','Fontsize', 15)
hLeg = legend({'$y_m = -1$','$y_m = +1$'},'Interpreter','latex','Fontsize', 15,'Location','north');
set(hLeg,'visible','off')
set(gca,'fontsize',15)
xlim([0 2])
ylim([0 1])

fig = figure; clf
box on
s2 = gscatter(DATA2(:,1)/pi,DATA2(:,2)/pi,DATA2(:,3));
s2(1).MarkerSize = 7; s2(1).Color = 'red';
s2(2).MarkerSize = 7; s2(2).Color = 'blue';
title('$\mbox{Markovian noise}$','Interpreter','latex','Fontsize', 15)
xlabel('$\theta_m/\pi$','Interpreter','latex','Fontsize', 15)
ylabel('$\phi_m/\pi$','Interpreter','latex','Fontsize', 15)
hLeg = legend({'$y_m = -1$','$y_m = +1$'},'Interpreter','latex','Fontsize', 15,'Location','north');
set(hLeg,'visible','off')
set(gca,'fontsize',15)
xlim([0 2])
ylim([0 1])

toc