clear all; close all;

%%
fhss;
time_freq_analysis(signal, -20, 100, 'Frequency Hopping Spread Spectrum');


%%
close all;

h = [1 zeros(1,16)];
SNR = 20;
modOrd = 4;
mu = 16;
N = 64;
ofdm_symbols = gen_ofdm(h, SNR, modOrd, mu, N, N);
time_freq_analysis(ofdm_symbols, 0, 2, 'OFDM');


