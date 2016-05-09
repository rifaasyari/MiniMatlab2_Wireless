clear all; close all;

%%

fhss_var = 0;
rngsettings = rng;

fhss;

pb = past_approach(freq_hopped_sig, -10, 100, 'FHSS (prev)');
time_freq_analysis(freq_hopped_sig, -10, 100, 'FHSS');

fhss_var = .25;
rng(rngsettings);
fhss;


[dummy, B1] = past_approach(freq_hopped_sig, -10, 100, 'FHSS noisy (prev)');
[dummy, B2] = time_freq_analysis(freq_hopped_sig, -10, 100, 'FHSS noisy');

p1 = sum(sum(B1 == pb)) / prod(size(pb))
p2 = sum(sum(B2 == pb)) / prod(size(pb))




%%
close all;

h = [1 zeros(1,16)];

modOrd = 4;
mu = 16;
N = 64;
rngsettings = rng;

SNR = 20;
ofdm_symbols = gen_ofdm(h, SNR, modOrd, mu, N, N);

pb = past_approach(ofdm_symbols, -10, 2, 'OFDM (prev)');
time_freq_analysis(ofdm_symbols, -10, 2, 'OFDM');


SNR = 6;
ofdm_symbols = gen_ofdm(h, SNR, modOrd, mu, N, N);
rng(rngsettings);
[dummy, B1] = past_approach(ofdm_symbols, -10, 2, 'OFDM noisy (prev)');
[dummy, B2] = time_freq_analysis(ofdm_symbols, -10, 2, 'OFDM noisy');

p1 = sum(sum(B1 == pb)) / prod(size(pb))
p2 = sum(sum(B2 == pb)) / prod(size(pb))
