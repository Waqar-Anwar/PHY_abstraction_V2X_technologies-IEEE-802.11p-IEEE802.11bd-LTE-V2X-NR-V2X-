% Currently these results are based on the PUSCH performmance with OFDM which uses almost simialr 
% PHY configration as defined by PSSCH. 

% This script is used to simulate the PER (Packet Error Ratio)

% Please Download first "Results_nr_awgn.mat" from GitHub and "channel_data.mat" from (https://www.dropbox.com/s/tx975gatia41mq1/channel_data.mat?dl=0)

b = load('Results_nr_awgn.mat'); % AWGN look-up table
chan = load('channel_data.mat'); % Channel gain recorded for 1000 channel relizations of Uran-crossing NLOS channel [M. Kahn, “V2V radio channel models,” IEEE 802.11-14/0259r0, Feb.2014] using PDSCH MATLAB 5G toolbox example [https://de.mathworks.com/help/5g/ug/nr-pdsch-throughput.html]

mcs = 21; % used MCS
ind =mcs+1;

beta = [2 2 2 2 2 2 2 2 2 2 10 10 10 10 10 10 10 42 42 42 42 42 42 42 42 42 42 42 42]; % Depending on the modulation 1 for BPSK, 2 for QPSK, 10 for 16-QAM, 42 for 64-QAM, 170 for 256-QAM obtained from SER (symbol error rate) experssions 

snr_d = 0:30;  % SNR at d considering the log normal path loss 


ICI = 0.000223;           % ICI calculated using eq. 31 of the paper

sigma_sym = abs(sinc(500*(1/60e3)))^2; % received signal power when the sampling is offset by Doppelr

for i=1:length(snr_d)
    
    pn = 10^(-snr_d(i)/10);
    ICI_plus_noise = (ICI + pn)/sigma_sym;
    
    snr_real = (abs(channel).^2) ./ICI_plus_noise;

        
    snr_awgn = b.snr_nr_awgn(ind,:);    % AWGN SINR table
    per_awgn = b.per_nr_awgn(ind,:);    % AWGN PER
    
    snr_ieesm=10*log10((beta(ind)/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta(ind))./sqrt(((2.*snr_real)./beta(ind)) + 1),[2 3]).^(-2))) - 1)); %Effective SINR mapping for used number of Resource blocks
    snr_ieesm(isinf(snr_ieesm)) = 100; %In case of high SINR matlab return +infinity and this replace it with some logical value
    loc_d = knnsearch(snr_awgn',snr_ieesm); % closest AWGN simulated points
    per_nr(i) = mean(per_awgn(loc_d));  %PER against AWGN SNR

end