function gen_q(Nt,M,K,H,Pt,n_channels,file_name_result)

eff_precoders = zeros(Nt,M,K,n_channels);
eff_rate = zeros(K,n_channels);
F = zeros(Nt,Nt,K);
W = zeros(Nr,Nr,K);
for nn = 1:n_channels
    nn
    Hn = H(:,:,:,nn);
    % Obtain F[k], W[k] with SVD

    q = zeros(Nt*M,K);
    rate = zeros(K,1);
    for k = 1:K
        k
        Hk = Hn(:,:,k);
        [U,S,V] = svd(Hk);
        p = waterfilling(Pt,1,diag(S));
        F(:,:,k) = U*sqrt(diag(p));
        W(:,:,k) = V;
        rate(k) = log2(det(eye(Nr) + pinv(W(:,:,k))*Hk*F(:,:,k)*F(:,:,k)'*Hk'*W(:,:,k)'));
    end

    eff_rate(:,nn) = rate;
end
save(file_name_result,'eff_precoders','eff_rate');