b = load('Results_nr_awgn.mat');
chan = load('estimation_mcs7_500.mat');
channel= chan.hest(1:1000,chan.psschIndices);

beta = [2 2 2 2 2 2 2 2 2 2 10 10 10 10 10 10 42 42 42 42];
data_rate_mcs = [1.40 1.87 2.31 2.94 3.60 4.42 5.11 6.07 6.94 8.10 8.10 8.83 9.71 10.79 12.14 13.88 16.19 19.43 24.25 32.38];
latency = [1.75 1.5 1.25 1 0.75 0.75 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.25 0.25 0.25 0.25 0.25 0.25 0.25];
nrb = [68 52 42 33 27 22 19 16 14 12 12 11 10 9 8 7 6 5 4 3];
valid_mcs = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 19 22 26];
mcs_len = 20;



d = 0.1:20:500;
snr_d = 23 - 47.86 -27.5*log10(d) - (-104) - 7;

% snr = 0:0.25:20;
clear data_rate latency_temp

ICI = 0.000223;

sigma_sym = abs(sinc(500*(1/15e3)))^2;

for i=1:length(snr_d)
    
    pn = 10^(-snr_d(i)/10);
    noise = (ICI + pn)/sigma_sym;
    %         noise = pn;
    
    %         snr_real = 10^(snr(i)/10).*(abs(channel).^2);
    
      estimat = reshape(channel,[length(channel(:,1)),length(channel(1,:))/10,10]);
      sigma_mmse(:,:)= mean(((abs(estimat).^2)./((abs(estimat).^2)+ noise )).^2,2);
      snr_real = 1./(1./(sqrt(sigma_mmse))-1);
    
    for ind =1:mcs_len
        
        snr_awgn = b.snr_lte_awgn(ind,:);
        per_awgn = b.per_lte_awgn(ind,:);
        
        snr_ieesm=10.*log10((beta(ind)/2).*(lambertw(exp(1).*((sum(exp(-snr_real'./beta(ind))./sqrt(((2.*snr_real')./beta(ind)) + 1))./10).^(-2))) - 1));
        snr_ieesm(isinf(snr_ieesm)) = 100;
        loc_d = knnsearch(snr_awgn',snr_ieesm');
        per = per_awgn(loc_d);
        data_rate(ind,:) = (1-per)*data_rate_mcs(ind);
        latency_temp(ind,:)=latency(ind)./((1-per)+0.001);
        
    end
    
    data_rate_nr(i) = mean(max(data_rate));
    latency_nr(i)= mean(min(latency_temp));
    
    %         snr_eesm = 10.*log10(-beta_esm*log(sum(exp(-snr_real'/beta_esm))./48));
    %         loc_d = knnsearch(snr_awgn',snr_eesm');
    %         per_11p_eesm_mcs2(i)=mean(per_awgn(loc_d));
    %
    %  for j=1:length(snr_real(:,1))
    %         loc1 = knnsearch(awgn_snr_mi', snr_real(j,:)'./beta_mi);
    %         pc_mi1 = mean(a.avg_mi(loc1));
    %         loc2 = knnsearch(a.avg_mi',pc_mi1);
    %         snr_mi(j) = 10.*log10(beta_mi.*awgn_snr_mi(loc2));
    %  end
    %
    %         loc_d = knnsearch(snr_awgn',snr_mi');
    %         per_11p_mi_mcs2(i)=mean(per_awgn(loc_d));
    
end



 
