clc;
clear all;
close all

%%  setup parameters
Nt = 8; % number of transmit antennas
Nr = 4;
Ns = Nr;
K = 64; % number of subcarriers

%% Load data
system = strcat('NrxNtxK=',num2str(Nr),'x',num2str(Nt),'x',num2str(K));
load(strcat(system,'.mat'));
system = strcat('NrxNtxK=',num2str(Nr),'x',num2str(Nt),'x',num2str(K));
C0 = Cbar(:,:,:,6);
%% ------------------------------------------------------------------------
SNR_dB = 10;
Pt = db2pow(SNR_dB);
n_channels = 5;

J_vec = [1,16:16:K];
T = length(theta);
markers = {'','+','s','d'}; lines = {'-',':'};

for rho = 0.25:0.25:0.75
    
    rate_sub = zeros(n_channels,length(J_vec));
    rate_all = zeros(n_channels,length(J_vec));
    MSE_sub = zeros(n_channels,length(J_vec));
    MSE_all = zeros(n_channels,length(J_vec));
    
    for jj = 1:length(J_vec)
        jj
        J = J_vec(jj);
        
        beam_all = zeros(T,K); beam_sub = zeros(T,J);
        for nn = 1:n_channels
            %nn
            Hn = H(:,:,:,nn);
            for k = 1:K
                Hk = Hn(:,:,k);
                [U,S,V] = svd(Hk);
                p = waterfilling(Pt,1,diag(S));
                F(:,:,k) = V(:,1:Ns)*sqrt(diag(p));
                W(:,:,k) = U(:,1:Ns);
                rate(k) = log2(real(det(eye(Ns) + 1/Ns * pinv(W(:,:,k))*Hk*F(:,:,k)*F(:,:,k)'*Hk'*W(:,:,k))));
            end
            
            [~, Omg] = mink(rate,J);
            [beam_ideal, beam_all_tmp, rate_all(nn,jj), MSE_all(nn,jj), beam_sub_tmp, rate_sub(nn,jj), MSE_sub(nn,jj)] = ...
                JCAS_design(Nt,Ns,K,Hn,C0,F,Pt,Omg,rho,at,T,Pd_theta);
            beam_all = beam_all + beam_all_tmp/n_channels;
            beam_sub = beam_sub + beam_sub_tmp/n_channels;
        end
        
        % average beams over all subcarriers
        beam_all_mean = mean(beam_all,2);
        beam_sub_mean = mean(beam_sub,2);
        
    end
    
    % Compute rate
    rate_mean_sub = mean(rate_sub,1);
    rate_mean_all = mean(rate_all,1);
    MSE_mean_sub = mean(MSE_sub,1);
    MSE_mean_all = mean(MSE_all,1);
    %% Plot gifures ========================================================================================
    system_Nrf = strcat('NtxNrxK=',num2str(Nt),'x',num2str(Nr),'x',num2str(K));
    
    %% plot rate vs J
    figure(3*(jj-1)+1)
    plot(J_vec, rate_mean_all, ':k','LineWidth',2,'MarkerSize',8); hold on;
    plot(J_vec, rate_mean_sub, '-r','LineWidth',2,'MarkerSize',8); hold on;
    xlabel('Number of JCAS subcarrier ($J$)','fontsize',12,'interpreter','latex');
    ylabel('Total achievable rate [bits/s/Hz]','fontsize',12,'interpreter','latex');
    xticks(J_vec)
    legend('All subcarrier',...
        'Subcarrier selection',...
        'Location','Best','fontsize',12,'interpreter','latex')
    title(strcat('Digital BF, ',system_Nrf));
    xlim([J_vec(1) J_vec(end)])
    
    
    %% plot rate vs J
    figure(3*(jj-1)+2)
    plot(J_vec, MSE_mean_all, ':k','LineWidth',2,'MarkerSize',8); hold on;
    plot(J_vec, MSE_mean_sub, '-r','LineWidth',2,'MarkerSize',8); hold on;
    xlabel('Number of JCAS subcarrier ($J$)','fontsize',12,'interpreter','latex');
    ylabel('Average MSE','fontsize',12,'interpreter','latex');
    legend('All subcarrier',...
        'Subcarrier selection',...
        'Location','Best','fontsize',12,'interpreter','latex')
    title(strcat('Digital BF, ',system_Nrf));
    xlim([J_vec(1) J_vec(end)])
    xticks(J_vec)
    
    %% plot rate vs MSE
    figure(3*(jj-1)+3)
    plot(MSE_mean_all, rate_mean_all, ':k','LineWidth',2,'MarkerSize',8); hold on;
    plot(MSE_mean_sub, rate_mean_sub, '-r','LineWidth',2,'MarkerSize',8); hold on;
    xlabel('Average MSE','fontsize',12,'interpreter','latex');
    ylabel('Achievable rate [bits/s/Hz]','fontsize',12,'interpreter','latex');
    legend('All subcarrier',...
        'Subcarrier selection',...
        'Location','Best','fontsize',12,'interpreter','latex')
    title(strcat('Digital BF, ',system_Nrf));
    %xlim([J_vec(1) J_vec(end)])
    %xticks(J_vec)
end
