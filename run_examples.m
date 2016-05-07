clear all; close all;

%%
fhss;
past_approach(signal, -20, 100, 'FHSS (previous work)');


%%
close all;

h = [1 zeros(1,16)];
SNR = 20;
modOrd = 4;
mu = 16;
N = 64;
ofdm_symbols = gen_ofdm(h, SNR, modOrd, mu, N, N);
past_approach(ofdm_symbols, 0, 2, 'OFDM (previous work)');


