# PHY abstraction of V2X technologies (IEEE-802.11p-IEEE802.11bd-LTE-V2X-NR-V2X)

This repository contains the code of PHY abstraction for current V2X technologies

Introduction:
PHY abstraction converts the fading SINR into AWGN equivalent SNR using effective SINR mapping techniques. After the conversion, effective SINR can be treated as AWGN SNR, which means the simulated performance in presence of the AWGN channel can be used to estimate performance under various fading conditions. This reasons for doing so are the following
1) To evaluate the performance of certain technology system-level simulations are used prior to its deployment for suitability analysis
2) The system-level simulation becomes infeasible if complete PHY is simulated per each node such as coding, constellation mapping, waveform generation, and passing waveform from the channel and reverse processing at the receiver. As these simulation require a huge amount of time and computations
3) To speed up simulation PHY abstraction is used, which map the received SINR to throughput or PER. 
4) As in the case of coded modulation no closed-form PER or data rate expression exists. Therefore, simulations are performed in presence of the AWGN channel for each MCS, and SINR vs. PER or Data rate tables are generated.
5) Using the PHY abstraction algorithms the fading SINR per symbol is converted to single AWGN equivalent SINR which then could be used to evaluate performance using the generated tablöes in step (4)

Code related details:
AWGN tables are generated using technology-specific MATLAB toolbox examples. Links are here
IEEE 802.11p: (https://de.mathworks.com/help/wlan/examples/802-11p-packet-error-rate-simulation-for-a-vehicular-channel.html)
IEEE 802.11bd: (https://de.mathworks.com/help/wlan/examples/802-11ac-packet-error-rate-simulation-for-8x8-tgac-channel.html)
LTE-V2X: (https://de.mathworks.com/help/lte/examples/release-14-v2x-sidelink-pssch-throughput.html)
NR-V2X: (https://de.mathworks.com/help/5g/ug/nr-pdsch-throughput.html)

Note: AWGN simulated tables are provided for technologies, which are required for PHY abstraction


Assumptions:
Ideal channel estimation, and perfect time and frequency synchronization is available. 


Channel Model:
V2X triger team defined Urbun Crossing NLOS channel model (M. Kahn, “V2V radio channel models” IEEE 802.11-14/0259r0, Feb.2014) is used, which are also provided in seprate MATLAB file.

