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

%% ------------------------------------------------------------------------
J_vec = [8,16,32];
plot_beampattern = 1; plot_rate = 1;

n_channels = 2;

SNR_vec = 0:2:12;
Pd_theta = Pd_theta;
T = length(theta);
markers = {'','+','s','d'}; lines = {'-',':'};

for jj = 1:length(J_vec)
    J = J_vec(jj);
    idx = 1;
    
    for rho = 0.5
        
        rate_sub = zeros(n_channels,length(SNR_vec));
        rate_all = zeros(n_channels,length(SNR_vec));
        MSE_sub = zeros(n_channels,length(SNR_vec));
        MSE_all = zeros(n_channels,length(SNR_vec));
        
        for ss = 1:length(SNR_vec)
            ss
            Pt = db2pow(SNR_vec(ss));
            C0 = Cbar(:,:,:,ss);
            
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
                    wk = W(:,1,k);
                    rate(k) = log2(real(det(eye(Ns) + 1/Ns * pinv(W(:,:,k))*Hk*F(:,:,k)*F(:,:,k)'*Hk'*W(:,:,k))));
                end
                
                [~, Omg] = mink(rate,J);
                [beam_ideal, beam_all_tmp, rate_all(nn,ss), MSE_all(nn,ss), beam_sub_tmp, rate_sub(nn,ss), MSE_sub(nn,ss)] = ...
                    JCAS_design(Nt,Ns,K,Hn,C0,F,Pt,Omg,rho,at,T,Pd_theta);
                beam_all = beam_all + beam_all_tmp/n_channels;
                beam_sub = beam_sub + beam_sub_tmp/n_channels;
            end
            
            % average beams over all subcarriers
            beam_all_mean = mean(beam_all,2);
            beam_sub_mean = mean(beam_sub,2);
            
        end
        
        %% plot beampattern
        if plot_beampattern == 1
            figure(3*(jj-1)+1)
            plot(theta*180/pi,Pd_theta,'--b','LineWidth',1);hold on;
            plot(theta*180/pi,mean(beam_all_tmp,2),':k','LineWidth',1);hold on;
            plot(theta*180/pi,mean(beam_sub_tmp,2),'-r','LineWidth',1);hold on;
            %plot(theta*180/pi,beam_all_mean,':k','LineWidth',2);hold on;
            %plot(theta*180/pi,beam_sub_mean,'-r','LineWidth',2);hold on;
            legend('Desired beampattern',...
                'All subcarrier',...
                'Subcarrier selection',...
                'Location','Best','fontsize',12,'interpreter','latex')
            xlim([-90, 90])
            xticks([-90:30:90])
            xlabel('Angles $(^{\circ})$ [dB]','fontsize',12,'interpreter','latex');
            ylabel('Normalized beampatter','fontsize',12,'interpreter','latex');
        end
        
        % Compute rate
        rate_mean_sub = mean(rate_sub,1);
        rate_mean_all = mean(rate_all,1);
        MSE_mean_sub = mean(MSE_sub,1);
        MSE_mean_all = mean(MSE_all,1);
        if plot_rate == 1
            %% Plot gifures ========================================================================================
            system_Nrf = strcat('NtxNrxK=',num2str(Nt),'x',num2str(Nr),'x',num2str(K));
            
            %% plot rate vs J
            figure(3*(jj-1)+2)
            plot(SNR_vec, rate_mean_all, strcat(':k',markers{idx}),'LineWidth',2,'MarkerSize',8); hold on;
            plot(SNR_vec, rate_mean_sub, strcat('-r',markers{idx}),'LineWidth',2,'MarkerSize',8); hold on;
            xlabel('SNR [dB]','fontsize',12,'interpreter','latex');
            ylabel('Total achievable rate [bits/s/Hz]','fontsize',12,'interpreter','latex');
            xticks(SNR_vec)
            legend('All subcarrier',...
                'Subcarrier selection',...
                'Location','Best','fontsize',12,'interpreter','latex')
            title(strcat('Digital BF, ',system_Nrf));
            xlim([SNR_vec(1) SNR_vec(end)])
            
            
            %% plot rate vs J
            figure(3*(jj-1)+3)
            plot(SNR_vec, MSE_mean_all, strcat(':k',markers{idx}),'LineWidth',2,'MarkerSize',8); hold on;
            plot(SNR_vec, MSE_mean_sub, strcat('-r',markers{idx}),'LineWidth',2,'MarkerSize',8); hold on;
            xlabel('SNR [dB]','fontsize',12,'interpreter','latex');
            ylabel('Average MSE','fontsize',12,'interpreter','latex');
            legend('All subcarrier',...
                'Subcarrier selection',...
                'Location','Best','fontsize',12,'interpreter','latex')
            title(strcat('Digital BF, ',system_Nrf));
            xlim([SNR_vec(1) SNR_vec(end)])
            xticks(SNR_vec)
        end
        
        idx = idx + 1;
    end
end