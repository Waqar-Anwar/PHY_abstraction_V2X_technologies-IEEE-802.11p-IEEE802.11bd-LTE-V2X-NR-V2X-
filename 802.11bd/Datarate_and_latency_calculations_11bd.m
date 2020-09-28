b = load('Results_11bd_awgn.mat');  % AWGN reference of 11bd cueves for first 10 MCS (0-9)
chan = load('channel_11bd.mat'); % fading channel
channel = chan.hest(1:1000,:,:); % first 1000 channel relaizations


P_tx = 23; % transmitted power
P_L0 = 47.86; % path loss at reference distance 
n_exp = 2.75; % path loss exponant
P_n = 104; % Gussion noise for 10 MHz bandwidth
NF = 9; % noise figure
G_tx = 3;% transmit antenna gain
G_rx = 3;% receive antenna gain
Pb = 300; % Number of bytes in a packet


n_dc_11bd = 52; % number of data carriers
r_11bd = [1/2 1/2 3/4 1/2 3/4 2/3 3/4 5/6 3/4 5/6]; % code rate of first 10 MCSs (0-9)
n_bps = [1 2 2 4 4 6 6 6 8 8]; %number of bits transmitted by each modulation schemes for first 10 MCSs (0-9)
t_ma = [8 8 8 8 8 4 4 4 4 4]; % midambel periodicity
t_pre_11bd = 80e-6; % preamble duratation
t_AIFS = 32e-6; % arbitrary inter-frame space
t_sym = 8e-6; % OFDM symbol duration



n_sym_11bd = ceil((Pb*8)./(n_dc_11bd.*r_11bd.*n_bps)); % number of data symbols 
n_ma = floor((n_sym_11bd-1)./t_ma); % number of midambles
t_tx_11bd = (t_pre_11bd + t_AIFS + t_sym.*n_sym_11bd + t_sym.*n_ma)*1e3; % transmission latency of each MCS in ms
Gamma_11bd = (Pb*8 ./t_tx_11bd)/1e3; % Data rates of each MCS in Mbps

beta = [1 2 2 10 10 42 42 42 170 170]; % beta values depending on the modulation order 1 for BPSK, 2 for QPSK, 10 for 16 QAM, 42 for 64 QAM, 170 for 256 QAM
% data_rate_mcs = [2.31 4.17 5.77 6.98 9.09 11.11 11.54 12.00 13.64 13.64];
% latency = [1.040 0.576 0.416 0.344 0.264 0.216 0.208 0.200 0.176 0.176];
mcs_len = 10; % number of MCSs



d = 0.1:20:520; % distance
snr_d = P_tx - P_L0 - n_exp*(10*log10(d)) - (-P_n) - NF + G_tx + G_rx;  % SNR at d considering the log normal path loss 

ICI = 3.2359e-05;           % ICI calculated using eq. 31 of the paper at carrier frequency of 5.9 GHz and 500 Hz Doppler

sigma_sym = abs(sinc(500*(1/156.25e3)))^2; % received signal power when the sampling time is offset by Doppelr
 
% snr = 0:0.25:20;
clear data_rate latency_temp

for i=1:length(snr_d)
    
    pn = 10^(-snr_d(i)/10);
    noise = (ICI + pn)/sigma_sym;
    
    snr_real = 1./((1./((abs(channel).^2)./((abs(channel).^2)+ noise)))-1); % received SINR per each data symbols and sub-carrier
    
    for ind =1:mcs_len
        
        snr_awgn = b.snr_11bd_awgn(ind,:);
        per_awgn = b.per_11bd_awgn(ind,:);
        
%        snr_ieesm=10.*log10((beta(ind)/2).*(lambertw(exp(1).*((sum(exp(-snr_real(:,:,1:n_sym_11bd(ind))'./beta(ind))./sqrt(((2.*snr_real')./beta(ind)) + 1))./52).^(-2))) - 1));
        snr_ieesm=10*log10((beta(ind)/2).*(lambertw(exp(1).*(mean(exp(-snr_real(:,:,1:n_sym_11bd(ind))./beta(ind))./sqrt(((2.*snr_real(:,:,1:n_sym_11bd(ind)))./beta(ind)) + 1),[2 3]).^(-2))) - 1)); %Effective SINR mapping
        snr_ieesm(isinf(snr_ieesm)) = 100;  % in case of positive infinity it will replace the value with 100
        loc_d = knnsearch(snr_awgn',snr_ieesm);
        per = per_awgn(loc_d);
        data_rate(ind,:) = (1-per)*Gamma_11bd(ind);  % effective data rate
        latency_temp(ind,:)=t_tx_11bd(ind)./((1-per)+0.0001); % transmission latency of selected mcs 0.0001 is used to avoid infinty in case of 1-per is 0
        
    end
    
    data_rate_11bd(i) = mean(max(data_rate)); % Select the MCS which deliver maximum data rates
    t_IAT_11bd(i)= mean(min(latency_temp)); % Select the MCS provide Minimum latency

    
end



 
