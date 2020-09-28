% Please first Download "Results_lte_awgn.mat" from GitHub and channel_data.mat from (https://www.dropbox.com/s/tx975gatia41mq1/channel_data.mat?dl=0)

b = load('Results_lte_awgn.mat'); % AWGN look-up table
chan = load('channel_data.mat'); % Channel gain recorded for 1000 channel relizations of Uran-crossing NLOS channel [M. Kahn, “V2V radio channel models,” IEEE 802.11-14/0259r0, Feb.2014] using PDSCH MATLAB 5G toolbox example [https://de.mathworks.com/help/5g/ug/nr-pdsch-throughput.html]
channel= chan.hest(1:1000,:,[1:4 6:7 9:10 12:13]); % Symbols used to transmit data


P_tx = 23; % transmitted power
P_L0 = 47.86; % path loss at reference distance 
n_exp = 2.75; % path loss exponant
P_n = 104; % Gussion noise for 10 MHz bandwidth
NF = 9; % noise figure
G_tx = 3;% transmit antenna gain
G_rx = 3;% receive antenna gain
Pb = 300; % Number of bytes in a packet

n_RB_lte = [86 66 54 41 34 28 23 20 18 16 14 14 12 11 10 9 8 8 7 7 6 6 6 5 5 5 4 4 4];  % number of required RBs to transmit a payload of 2400 bits Table 7.1.7.2.1-1 of TS 36.213
t_fr = 10e-3; % frame duration
t_sub_fr = 1e-3; % sub-frame duration
n_RB_fr_lte = 408;  % available RBs in a frame for data communication Tab.2 in [L. Tung and M. Gerla, “LTE resource scheduling for vehicular safety applications,”in 10th Annual Conference on Wireless On-demand Network Systems and Services (WONS), Mar. 2013, pp. 116–118.]


Gamma_lte = ((Pb*8*n_RB_fr_lte)./(n_RB_lte.*t_fr))./1e6; % Mbps
t_tx_lte = ceil(n_RB_lte./(n_RB_fr_lte/10)).*(t_sub_fr*1e3); % transmission latency in ms

beta = [2 2 2 2 2 2 2 2 2 2 2 10 10 10 10 10 10 10 10 42 42 42 42 42 42 42 42 42 42]; % Depending on the modulation 1 for BPSK, 2 for QPSK, 10 for 16-QAM, 42 for 64-QAM obtained from SER (symbol error rate) experssions 
mcs_len = 10; % number of MCSs



d = 0.1:20:520; % distance
snr_d = P_tx - P_L0 - n_exp*(10*log10(d)) - (-P_n) - NF + G_tx + G_rx;  % SNR at d considering the log normal path loss 

ICI = 3.6e-03;           % ICI calculated using eq. 31 of the paper at carrier frequency of 5.9 GHz and 500 Hz Doppler

sigma_sym = abs(sinc(500*(1/15e3)))^2; % received signal power when the sampling time is offset by Doppelr
 
% snr = 0:0.25:20;
clear data_rate latency_temp

for i=1:length(snr_d)
    
    pn = 10^(-snr_d(i)/10);
    ICI_plus_noise = (ICI + pn)/sigma_sym;
    
  
    for ind =1:mcs_len
        
        snr_awgn = b.snr_lte_awgn(ind,:);
        per_awgn = b.per_lte_awgn(ind,:);
        
        sigma_mmse = mean(((abs(channel(:,1:12*n_RB_lte(ind),:)).^2)./((abs(channel(:,1:12*n_RB_lte(ind),:)).^2)+ ICI_plus_noise )).^2,2);
        snr_real = 1./(1./(sqrt(sigma_mmse))-1); % Received SINR per each DFT-s-OFDM symbol 
        snr_ieesm=10*log10((beta(ind)/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta(ind))./sqrt(((2.*snr_real)./beta(ind)) + 1),3).^(-2))) - 1)); %Effective SINR mapping          
        snr_ieesm(isinf(snr_ieesm)) = 100;  % in case of positive infinity it will replace the value with 100

        loc_d = knnsearch(snr_awgn',snr_ieesm);
        per = per_awgn(loc_d);
        data_rate(ind,:) = (1-per)*Gamma_lte(ind);  % effective data rate
        latency_temp(ind,:)=t_tx_lte(ind)./((1-per)+0.0001); % transmission latency of selected mcs 0.0001 is used to avoid infinty in case of 1-per is 0
        
    end
    
    data_rate_lte(i) = mean(max(data_rate)); % Select the MCS which deliver maximum data rates
    t_IAT_lte(i)= mean(min(latency_temp)); % Select the MCS provide Minimum latency

    
end



 
