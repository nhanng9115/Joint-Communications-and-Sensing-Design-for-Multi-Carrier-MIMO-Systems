function [q0, rate] = optimize_q(Nt,M,K,h,x0,q0,SNR)




%opt_solver = sdpsettings('solver','mosek','verbose',0);

convergence = 0;
obj_old = 0;
count = 0;
count_max = 100;
% while convergence == 0
    
    cvx_solver sedumi
    cvx_begin quiet
    
    %% variables
    variable q(Nt*M,K)
    variable x
    variable t(M,K)
    %q = sdpvar(Nt*M,K,'full','complex');
    %x = sdpvar(M,K,'full','real');
    
    %% objective function
    obj = 0;
    for k = 1:K
        for m = 1:M
            obj = obj + t(m,k);
        end
    end
    maximize obj
    
    %% constraints
    subject to

    
    % second constrant
    for k = 1:K
        norm(q(:,k))^2 <= SNR; % normalize to get SNR as the max power
        for m = 1:M
            hmk = h(:,m,k);
            Hmk = hmk*hmk';
            
            % construct Hi_hat and Hi_bar
            Hmk_hat = zeros(M*Nt,M*Nt);
            Hmk_hat(Nt*(m-1)+1:Nt*m,Nt*(m-1)+1:Nt*m) = Hmk;
            
            Hmk_bar = kron(eye(M),Hmk);
            Hmk_bar(Nt*(m-1)+1:Nt*m,Nt*(m-1)+1:Nt*m) = zeros(Nt,Nt);
            
            F_qol = q0(:,k)'*Hmk_hat*q0(:,k)/x0(m,k)^2 * x(m,k) - 2*real(q0(:,k)'*Hmk_hat*q(:,k))/x0(m,k); % (17)
            q(:,k)'*Hmk_bar*q(:,k) + 1 + F_qol <= 0; % (18)
            x(m,k) <= 2^t - 1;
        end
    end
    
    cvx_end
    %% Solving
    %optimize(F, -obj, opt_solver);
    
    %% update variables
    q0 = double(q);
    x0 = double(x);
    
    %% compute objective value
    obj_new = double(obj);
    if abs(obj_new - obj_old) <= 1e-3 || count >= count_max
        convergence = 1;
    else
        obj_old = obj_new;
        count = count + 1;
    end
% end
count
rate_all = log2(1 + x0);
rate = sum(rate_all,1); % rate for subcarriers
end