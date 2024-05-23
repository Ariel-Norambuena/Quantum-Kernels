tic

%% Physical Parameters
% You can change this part to generate different data sets
% 
J = 1;
h = 0.6*J;
alpha = 0.5;

%% Generation of data
Nspins = 2;
DATA = GenerateDataDQPT(J,h,alpha,Nspins(k));
string1 = 'DATA_N=';
string2 = num2str(Nspins);
string3 = '_theta_phi_DQPT.txt';
name = strcat(string1,string2,string3);
writematrix(DATA,name);

toc