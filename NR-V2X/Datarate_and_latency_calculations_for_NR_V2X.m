% Currently these results are based on the PUSCH performmance with OFDM which uses almost simialr 
% PHY configration as defined by PSSCH. 

% Please Download first Results_nr_awgn.mat from GitHub and channel_data.mat from (https://www.dropbox.com/s/tx975gatia41mq1/channel_data.mat?dl=0)

b = load('Results_nr_awgn.mat'); % AWGN look-up table
chan = load('channel_data.mat'); % Channel gain recorded for 1000 channel relizations of Uran-crossing NLOS channel [M. Kahn, �V2V radio channel models,� IEEE 802.11-14/0259r0, Feb.2014] using PDSCH MATLAB 5G toolbox example [https://de.mathworks.com/help/5g/ug/nr-pdsch-throughput.html]
channel= chan.hest(1:1000,:,2:13); % Symbols used to transmit data

P_tx = 23; % transmitted power
P_L0 = 47.86; % path loss at reference distance 
n_exp = 2.75; % path loss exponant
P_n = 104; % Gussion noise for 10 MHz bandwidth
NF = 9; % noise figure
G_tx = 3;% transmit antenna gain
G_rx = 3;% receive antenna gain
Pb = 300; % number of bytes in a packet
n_RE_NR = 105; % 3 DMRS symbols multiplexed with PSSCH have 18 RE (resource element) + 2 PSCCH multiplexed with PSSCH have 12 RE + 1 DMRS symbol multiplexed with PSCCH and PSSCH will have 3 RE and + 6 PSSCH symbole will have 72 RE; 3GPP 21.916, Fig.1
n_RB_slot_NR= 11; % 60 KHz spacing and 10 MHz bandwidth, TS 38.101-1   Table 5.3.2-1
t_slot =0.00025;  % for 60 KHz spacing, 38.101-1 TS 38.211  Table 4.3.2-1
n_bpm =[0.2344 0.3066 0.3770 0.4902 0.6016 0.7402 0.8770 1.0273 1.1758 1.3262 1.3281 1.4766 1.6953 1.9141 2.1602 2.4063 2.5703 2.5664 2.7305 3.0293 3.3223 3.6094 3.9023 4.2129 4.5234 4.8164 5.1152 5.3320 5.5547]; % from 3GPP TS 38.214 V16.2.0, Table 5.1.3.1-1

n_RB_NR = ceil((Pb*8)./(n_RE_NR.*n_bpm));
Gamma_NR = ((Pb*8*n_RB_slot_NR)./(n_RB_NR.*t_slot))./1e6; % Mbps
t_tx_NR = ceil(n_RB_NR/n_RB_slot_NR).*(t_slot*1e3); % transmission latency in ms

mcs_len = length(n_bpm) -1; % first 27 MCS are considered

beta = [2 2 2 2 2 2 2 2 2 2 10 10 10 10 10 10 10 42 42 42 42 42 42 42 42 42 42 42 42]; % Depending on the modulation 1 for BPSK, 2 for QPSK, 10 for 16-QAM, 42 for 64-QAM, 170 for 256-QAM obtained from SER (symbol error rate) experssions 

d = 0.1:20:520; % distance
snr_d = P_tx - P_L0 - n_exp*(10*log10(d)) - (-P_n) - NF + G_tx + G_rx;  % SNR at d considering the log normal path loss 

% snr = 0:0.25:20;
clear data_rate latency_temp

ICI = 0.000223;           % ICI calculated using eq. 31 of the paper at carrier frequency of 5.9 GHz and 500 Hz Doppler

sigma_sym = abs(sinc(500*(1/60e3)))^2; % received signal power when the sampling is offset by Doppelr

for i=1:length(snr_d)
    
    pn = 10^(-snr_d(i)/10);
    ICI_plus_noise = (ICI + pn)/sigma_sym;
    
    snr_real = (abs(channel).^2) ./ICI_plus_noise;
    
    
    for ind =1:mcs_len
        
        snr_awgn = b.snr_nr_awgn(ind,:);    % AWGN SINR table
        per_awgn = b.per_nr_awgn(ind,:);    % AWGN PER
        
        snr_ieesm=10*log10((beta(ind)/2).*(lambertw(exp(1).*(mean(exp(-snr_real(:,1:(12*n_RB_NR(ind)),:)./beta(ind))./sqrt(((2.*snr_real(:,1:(12*n_RB_NR(ind)),:))./beta(ind)) + 1),[2 3]).^(-2))) - 1)); %Effective SINR mapping
        snr_ieesm(isinf(snr_ieesm)) = 100; %In case of high SINR matlab return +infinity and this replace it with some logical value
        loc_d = knnsearch(snr_awgn',snr_ieesm); % closest AWGN simulated points
        per = per_awgn(loc_d);  %PER against AWGN SNR
        data_rate(ind,:) = (1-per)*Gamma_NR(ind);  % data rate calcuation for the considered MCS
        latency_temp(ind,:)= t_tx_NR(ind)./((1-per)+ 0.0001); % data rate calcuation for the considered MCS and 0.0001 is just used to avoid infinity
        
    end
    
    data_rate_nr(i) = mean(max(data_rate));  % Select the MCS which deliver maximum data rates 
    t_IAT_nr(i)= mean(min(latency_temp));  % Select the MCS provide Minimum latency
end



 
