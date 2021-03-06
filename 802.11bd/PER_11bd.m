% This script is used to simulate the PER (Packet Error Ratio) for 802.11bd

% Please Download first "Results_nr_awgn.mat" from GitHub and "channel_data.mat" from (https://www.dropbox.com/s/tx975gatia41mq1/channel_data.mat?dl=0)

% Please Download first Results_11bd_awgn.mat from GitHub and channel_11bd.mat from (which can be downloaded here: https://www.dropbox.com/s/3522ashf40qjz1a/channel_11bd.mat?dl=0)

b = load('Results_11bd_awgn.mat');  % AWGN reference of 11bd cueves for first 10 MCS (0-9)
chan = load('channel_11bd.mat'); % fading channel
channel = chan.hest(1:1000,:,:); % first 1000 channel relaizations

mcs = 5; % used MCS 
ind = mcs+1;

beta = [1 2 2 10 10 42 42 42 170 170]; % beta values depending on the modulation order 1 for BPSK, 2 for QPSK, 10 for 16 QAM, 42 for 64 QAM, 170 for 256 QAM

snr_d = 0:30;  % SNR at d considering the log normal path loss 

ICI = 3.2359e-05;           % ICI calculated using eq. 31 of the paper at carrier frequency of 5.9 GHz and 500 Hz Doppler

sigma_sym = abs(sinc(500*(1/156.25e3)))^2; % received signal power when the sampling time is offset by Doppelr

for i=1:length(snr_d)
    
    pn = 10^(-snr_d(i)/10);
    ICI_plus_noise = (ICI + pn)/sigma_sym;
    
    snr_real = (abs(channel).^2) ./ICI_plus_noise;

        
    snr_awgn = b.snr_11bd_awgn(ind,:);    % AWGN SINR table
    per_awgn = b.per_11bd_awgn(ind,:);    % AWGN PER
    
    snr_ieesm=10*log10((beta(ind)/2).*(lambertw(exp(1).*(mean(exp(-snr_real./beta(ind))./sqrt(((2.*snr_real)./beta(ind)) + 1),[2 3]).^(-2))) - 1)); %Effective SINR mapping for used number of symbols
    snr_ieesm(isinf(snr_ieesm)) = 100; %In case of high SINR matlab return +infinity and this replace it with some logical value
    loc_d = knnsearch(snr_awgn',snr_ieesm); % closest AWGN simulated points
    per_11bd(i) = mean(per_awgn(loc_d));  %PER against AWGN SNR

end