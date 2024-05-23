function rho = SolvesMasterEquation(L,rho_0,tspan)

N = size(L,1);
TOL = 1e-6*ones(1,N);
options = odeset('RelTol',1e-6,'AbsTol',TOL);
[~,Y] = ode45(@fun,tspan,rho_0,options);
SOL = Y;
% Reshape solution
dim = sqrt(N);
rho = cell(1,length(tspan));
for i=1:length(tspan)
    rho{i} = reshape(SOL(i,:),dim,dim);
end

    function dy = fun(~,y)
        dy = L*y;
    end

end