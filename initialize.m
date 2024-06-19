function [q0, x0] = initialize(M,Nt,K,h,SNR)
%% Random initial point ----------------------------------------
norm_min = 100;
for i_rand = 1:1
    q0_tmp = sqrt(SNR/2) * (randn(M*Nt,K)+ 1j*randn(M*Nt,K));
    if norm(q0_tmp) <= norm_min
        q0 = q0_tmp;
        norm_min = norm(q0_tmp);
    end
end

x0 = zeros(M,K);
for k = 1:K
    for m = 1:M
        hmk = h(:,m,k);
        Hmk = hmk*hmk';
        
        % construct Hi_hat and Hi_bar
        Hmk_hat = zeros(M*Nt,M*Nt);
        Hmk_hat(Nt*(m-1)+1:Nt*m,Nt*(m-1)+1:Nt*m) = Hmk;
        
        Hmk_bar = kron(eye(M),Hmk);
        Hmk_bar(Nt*(m-1)+1:Nt*m,Nt*(m-1)+1:Nt*m) = zeros(Nt,Nt);
        
        x0(m,k) = real(q0(:,k)'*Hmk_hat*q0(:,k) / (q0(:,k)'*Hmk_bar*q0(:,k) + 1));
    end
end
end