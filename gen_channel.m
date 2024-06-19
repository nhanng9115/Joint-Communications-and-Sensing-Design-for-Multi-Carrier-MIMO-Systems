clear all;
clc;

Nt = 4; % number of transmit antennas
Nr = 2;
K = 2;
Pt = db2pow(10);
channel_type = 0;

f0 = 2e6; DeltaF = 100e3; fmax = f0 + (K-1)*DeltaF;
system_config = strcat('NrxNtxK=',num2str(Nr),'x',num2str(Nt),'x',num2str(K));

n_trial = 1;
H = zeros(Nr,Nt,K,n_trial);
theta = pi.*rand(n_trial,1) - pi/2;
phi = pi.*rand(n_trial,1) - pi/2;
for ii = 1:n_trial
    for k = 1:K
        fk = f0 + (k-1)*DeltaF;
        if channel_type == 0
            H(:,:,k,ii) = 1/sqrt(2) * (randn(Nr,Nt) + 1i*randn(Nr,Nt));
        else
            at = exp( -1i * pi * fk/fmax * sin(theta(ii)) .* [0:Nt-1]' );
            ar = exp( -1i * pi * fk/fmax * sin(phi(ii)) .* [0:Nr-1]' );
            H(:,:,k,ii) = ar*at';
        end
    end
end

%% Radar detection angles
delta = pi/180;
theta = -pi/2:delta:pi/2;
target_DoA = [-pi/3,-pi/6,pi/6,pi/3];
% target_DoA = [-pi/3,0,pi/3];

%% Ideal beampattern design
beam_width = 9;
l = ceil((target_DoA+pi/2*ones(1,length(target_DoA)))/(delta)+ones(1,length(target_DoA)));
Pd_theta = zeros(length(theta),1);
for ii = 1:length(target_DoA)
    Pd_theta(l(ii)-(beam_width-1)/2:l(ii)+(beam_width-1)/2,1) = ones(beam_width,1);
end

SNR_vec = 0:2:12;
at = zeros(Nt,length(theta),K);
Cbar = zeros(Nt,Nt,K,length(SNR_vec));
for k = 1:K
    fk = f0 + (k-1)*DeltaF;
    for tt = 1:length(theta)
        at(:,tt,k) = exp( -1i * pi * fk/fmax * sin(theta(tt)) .* [0:Nt-1]' );
    end
    for ss = 1:length(SNR_vec)
        Pt = db2pow(SNR_vec(ss));
        Cbar(:,:,k,ss) = solve_Cbar(Pd_theta,Nt,at(:,:,k),theta,Pt);
    end
end

save(strcat(system_config,'.mat'));