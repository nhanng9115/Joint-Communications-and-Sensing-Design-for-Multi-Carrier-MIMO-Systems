function Q = opt_radcom(Nt,Ns,rho,Cbar,X1,Pt)

rho1 = (1-rho);

% Create the problem structure.
MO = obliquecomplexfactory(Nt,Ns);
problem.M = MO;

% Define the problem cost function and its gradient.

problem.cost = @cost;
    function f = cost(X)
        f = rho*norm(X*X' - Cbar,'fro')^2 + rho1*norm(X - X1,'fro')^2;
    end
problem.grad = @(X) problem.M.egrad2rgrad(X,egrad(X));
    function g = egrad(X)
        g = 4*rho*(X*X' - Cbar)*X + 2*rho1*(X - X1);
    end

% checkgradient(problem);
% warning('off', 'manopt:getHessian:approx');

% Execute the optimization
Qtmp = conjugategradient(problem);
Q = sqrt(Pt/Ns)*Qtmp;

% check power constraint
if abs(norm(Q,'fro')^2 - Pt) > 1e-4
    error('check power constraint!!!!!!!!!');
end

%test
% term1 = rho*norm(Q*Q' - Cbar,'fro')^2;
% term2 = rho1*norm(Q - X1,'fro')^2;
% ratio = term1 / term2;

end