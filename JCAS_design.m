function [beam_ideal, beam_all, rate_all, MSE_all, beam_sub, rate_sub, MSE_sub] = JCAS_design(Nt,Ns,K,H,Cbar,F0,Pt,Omg,rho,at,T,Pd_theta)
J = length(Omg);
Nr = size(H,1);
%% digital
F_all = zeros(Nt,Ns,K);
F_sub = F0;

beam_ideal = zeros(T,K);
beam_all = zeros(T,K);
beam_sub = zeros(T,J);
% figure

j = 1;
for k = 1:K
    beam_ideal(:,k) = real(diag(at(:,:,k)'*Cbar(:,:,k)*at(:,:,k))/Pt);
    Fbar = F0(:,:,k);
    F_all(:,:,k) = opt_radcom(Nt,Ns,rho,Cbar(:,:,k),Fbar,Pt);
    beam_all(:,k) = real(diag(at(:,:,k)'*F_all(:,:,k)*F_all(:,:,k)'*at(:,:,k))/Pt);
        
    % subcarrier selection: only use opt_radcom for k in Omg
    if ismember(k,Omg)
        F_sub(:,:,k) = F_all(:,:,k);
        beam_sub(:,j) = real(diag(at(:,:,k)'*F_sub(:,:,k)*F_sub(:,:,k)'*at(:,:,k))/Pt);
        j = j + 1;
    end
    
    
end
    
%% compute rate
rate_all = 0; rate_sub = 0;
for k = 1:K
    Hk = H(:,:,k);
    [Wk_all,~,~] = svd(Hk*F_all(:,:,k));
    [Wk_sub,~,~] = svd(Hk*F_sub(:,:,k));
    rate_all = rate_all + 1/K*log2(real(det(eye(Ns) + 1/Ns * pinv(Wk_all)*Hk*F_all(:,:,k)*F_all(:,:,k)'*Hk'*Wk_all)));
    rate_sub = rate_sub + 1/K*log2(real(det(eye(Ns) + 1/Ns * pinv(Wk_sub)*Hk*F_sub(:,:,k)*F_sub(:,:,k)'*Hk'*Wk_sub)));
end

%% Compute MSE
MSE_sub = norm(Pd_theta - mean(beam_sub,2))^2/T;
MSE_all = norm(Pd_theta - mean(beam_all,2))^2/T;

end % EOF