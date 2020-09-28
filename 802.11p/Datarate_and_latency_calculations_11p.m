% This script could be used to simulate Data rate and transmission latency
% for 802.11p

% Please Download first Results_11p_awgn.mat from GitHub and channel_11p.mat from (which can be downloaded here: https://www.dropbox.com/s/98lxwzqvzzoxp9l/channel_11p.mat?dl=0)

b = load('Results_11p_awgn.mat');  % AWGN reference cueves of 11p
chan = load('channel_11p.mat'); % fading channel
channel = chan.hest(1:1000,:,:); % first 1000 channel relaizations


P_tx = 23; % transmitted power
P_L0 = 47.86; % path loss at reference distance 
n_exp = 2.75; % path loss exponant
P_n = 104; % Gussion noise for 10 MHz bandwidth
NF = 9; % noise figure
G_tx = 3;% transmit antenna gain
G_rx = 3;% receive antenna gain
Pb = 300; % Number of bytes in a packet


n_dc_11p = 48; % number of data carriers
r_11p = [1/2 3/4 1/2 3/4 1/2 3/4 2/3 3/4]; % code rate of 11p MCSs
n_bps = [1 1 2 2 4 4 6 6]; %number of bits transmitted by each modulation schemes of 11p MCSs
t_pre_11p = 40e-6; % preamble duratation
t_AIFS = 32e-6; % arbitrary inter-frame space
t_sym = 8e-6; % OFDM symbol duration
beta = [1 1 2 2 10 10 42 42]; % beta values depending on the modulation order 1 for BPSK, 2 for QPSK, 10 for 16 QAM, 42 for 64 QAM



n_sym_11p = ceil((Pb*8)./(n_dc_11p.*r_11p.*n_bps)); % number of data symbols 
t_tx_11p = (t_pre_11p + t_AIFS + t_sym.*n_sym_11p)*1e3; % transmission latency of each MCS in ms
Gamma_11p = (Pb*8 ./t_tx_11p)/1e3; % Data rates of each MCS in Mbps

mcs_len = 8; % number of MCSs



d = 0.1:20:520; % distance
snr_d = P_tx - P_L0 - n_exp*(10*log10(d)) - (-P_n) - NF + G_tx + G_rx;  % SNR at d considering the log normal path loss 

ICI = 3.2359e-05;           % ICI calculated using eq. 31 of the paper at carrier frequency of 5.9 GHz and 500 Hz Doppler

sigma_sym = abs(sinc(500*(1/156.25e3)))^2; % received signal power when the sampling time is offset by Doppelr
 
% snr = 0:0.25:20;
clear data_rate latency_temp

for i=1:length(snr_d)
    
    pn = 10^(-snr_d(i)/10);
    ICI_plus_noise = (ICI + pn)/sigma_sym;
    
    snr_real = (abs(channel).^2) ./ICI_plus_noise; % received SINR per sub-carrier of each OFDM symbol
    
    for ind =1:mcs_len
        
        snr_awgn = b.snr_11p_awgn(ind,:);
        per_awgn = b.per_11p_awgn(ind,:);
        
        snr_ieesm=10*log10((beta(ind)/2).*(lambertw(exp(1).*(mean(exp(-snr_real(:,:,1:n_sym_11p(ind))./beta(ind))./sqrt(((2.*snr_real(:,:,1:n_sym_11p(ind)))./beta(ind)) + 1),[2 3]).^(-2))) - 1)); %Effective SINR mapping
        snr_ieesm(isinf(snr_ieesm)) = 100;  % in case of positive infinity it will replace the value with 100
        loc_d = knnsearch(snr_awgn',snr_ieesm);
        per = per_awgn(loc_d);
        data_rate(ind,:) = (1-per)*Gamma_11p(ind);  % effective data rate
        latency_temp(ind,:)=t_tx_11p(ind)./((1-per)+0.0001); % transmission latency of selected mcs 0.0001 is used to avoid infinty in case of 1-per is 0
        
    end
    
    data_rate_11p(i) = mean(max(data_rate)); % Select the MCS which deliver maximum data rates
    t_IAT_11p(i)= mean(min(latency_temp)); % Select the MCS provide Minimum latency
 
end



 
