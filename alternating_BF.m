function [Qall,Qsub,V,Omg] = alternating_BF(Nt,Nr,M,K,Nrf,H,Pt,Cbar,rho,J)
% design effective precoders (for both fully digital and HBF) and digital
% combiner
sigma2 = 1;
sumrate0 = 0; rate_list = [];
converge = 0;
Q0 = zeros(Nt,M,K);
eps = 1e-2;

iter = 0;
while converge == 0
    
    % Obtain effective precoder by solving sub-1
    if iter == 0
        % Initialization
        Q0 = sqrt(1/4) * (randn(Nt,M,K)+ 1j*randn(Nt,M,K));
    else
        Q0 = opt_precoder(Nt,Nr,M,K,H,Pt,q0,x0,V);
    end
    
    % Solving digital combiner
    for k = 1:K
        for m = 1:M
            gmmk = H(:,:,m,k)*Q0(:,m,k);
            Gmk = [];
            for j = 1:M
                gmjk = H(:,:,m,k)*Q0(:,j,k);
                Gmk = [Gmk,gmjk];
            end
            V(:,m,k) = (Gmk*Gmk' + sigma2*eye(Nr))^(-1) * gmmk;
        end
    end
    % compute rate
    [sumrate, SINR, rate] = compute_rate(H,Nt,M,K,zeros(Nt,Nrf),Q0,V);
    sumrate
    rate_list = [rate_list, sumrate];
    if abs(sumrate - sumrate0) <= eps
        converge = 1;
    else
        sumrate0 = sumrate;
        q0 = reshape(Q0,[Nt*M,K]);
        x0 = SINR;
        iter = iter + 1;
    end
end
figure
plot([1:length(rate_list)],rate_list);

% Obtain JCAS precoder by solving sub-2
Qall = zeros(Nt,M,K);
Qsub = zeros(Nt,M,K);
rate_K = sum(rate,1);
[~, Omg] = mink(rate_K,J);
for k = 1:K
    Qall(:,:,k) = opt_radcom(Nt,M,rho,Cbar(:,:,k),Q0(:,:,k),Pt);
    if ismember(k,Omg)
        Qsub(:,:,k) = Qall(:,:,k);
    else
        Qsub(:,:,k) = Q0(:,:,k);
    end
end


end % EOF