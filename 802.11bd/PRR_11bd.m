% This script could be used to simulate the PRR (Packet Reception Ratio)
% for 802.11bd

% Please Download first Results_11bd_awgn.mat from GitHub and channel_11bd.mat from (which can be downloaded here: https://www.dropbox.com/s/3522ashf40qjz1a/channel_11bd.mat?dl=0)

b = load('Results_11bd_awgn.mat');  % AWGN reference of 11bd cueves for first 10 MCS (0-9)
chan = load('channel_11bd.mat'); % fading channel
channel = chan.hest(1:1000,:,:); % first 1000 channel relaizations

mcs = 1; % used MCS 
ind = mcs+1;

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
beta = [1 2 2 10 10 42 42 42 170 170]; % beta values depending on the modulation order 1 for BPSK, 2 for QPSK, 10 for 16 QAM, 42 for 64 QAM, 170 for 256 QAM

d = 0.1:20:520; % distance
snr_d = P_tx - P_L0 - n_exp*(10*log10(d)) - (-P_n) - NF + G_tx + G_rx;  % SNR at d considering the log normal path loss 

ICI = 3.2359e-05;           % ICI calculated using eq. 31 of the paper at carrier frequency of 5.9 GHz and 500 Hz Doppler

sigma_sym = abs(sinc(500*(1/156.25e3)))^2; % received signal power when the sampling time is offset by Doppelr
 
% snr = 0:0.25:20;


for i=1:length(snr_d)
    
    pn = 10^(-snr_d(i)/10);
    ICI_plus_noise = (ICI + pn)/sigma_sym;
    
    snr_real = (abs(channel).^2) ./ICI_plus_noise;

        
    snr_awgn = b.snr_11bd_awgn(ind,:);    % AWGN SINR table
    per_awgn = b.per_11bd_awgn(ind,:);    % AWGN PER
    
    snr_ieesm=10*log10((beta(ind)/2).*(lambertw(exp(1).*(mean(exp(-snr_real(:,:,1:n_sym_11bd(ind))./beta(ind))./sqrt(((2.*snr_real(:,:,1:n_sym_11bd(ind)))./beta(ind)) + 1),[2 3]).^(-2))) - 1)); %Effective SINR mapping
    snr_ieesm(isinf(snr_ieesm)) = 100;  % in case of positive infinity it will replace the value with 100
    loc_d = knnsearch(snr_awgn',snr_ieesm); % closest AWGN simulated points
    prr_11bd(i) = 1-mean(per_awgn(loc_d));  %PER against AWGN SNR

end