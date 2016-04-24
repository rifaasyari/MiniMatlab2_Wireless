% Frequency Hopping Spread Spectrum
clc;
clear all;
close all;

% Generation of bits
N=20; 
s=randint(1,N);   % Generating N bits
signal=[];  
carrier=[];
T=180;
t=0:2*pi/(T-1):2*pi;     % Creating T samples for one cosine 
subplot(5,1,1);
stem(s);
set(gca,'XTick',1:N,'XLim',[0.66 N+0.33]);    %Setting axis limits and scale for the graph
title('Binary information');
ylabel('Amplitude(V)');
xlabel('Time(s)');

% Generation of bit pattern/bit stream(Polar NRZ type)
for k=1:N
     if s(1,k)==0
         sig=-ones(1,T);    % T no.of MINUS ONES for bit 0
     else
         sig=ones(1,T);     % T no.of ONES for bit 1
     end
     c=cos(t);   
     carrier=[carrier c];
     signal=[signal sig];
end
subplot(5,1,2);
stairs(signal);
set(gca,'XTick',0:T:N*T,'XLim',[-60 N*T+60],'YLim',[-1.5 1.5]);      %Setting axis limits and scale for the graph
title('Original Bit Sequence');
ylabel('Amplitude(V)');
xlabel('Time(s)');
 
% BPSK Modulation of the signal
bpsk_sig=signal.*carrier;   % Modulating the signal
subplot(5,1,3);
plot(bpsk_sig);
set(gca,'XTick',0:T:N*T,'XLim',[-60 N*T+60],'YLim',[-1.5 1.5]);
title('BPSK Modulated Signal');
ylabel('Amplitude(V)');
xlabel('Time(s)');
 
% Preparation of 6 new carrier frequencies
t1=0:2*pi/8:2*pi;
t2=0:2*pi/9:2*pi;
t3=0:2*pi/17:2*pi;
t4=0:2*pi/35:2*pi;
t5=0:2*pi/89:2*pi;
t6=0:2*pi/179:2*pi;
c1=cos(t1);
c1=[c1 c1 c1 c1 c1 c1 c1 c1 c1 c1 c1 c1 c1 c1 c1 c1 c1 c1 c1 c1];
c2=cos(t2);
c2=[c2 c2 c2 c2 c2 c2 c2 c2 c2 c2 c2 c2 c2 c2 c2 c2 c2 c2] ;
c3=cos(t3);
c3=[c3 c3 c3 c3 c3 c3 c3 c3 c3 c3] ;
c4=cos(t4);
c4=[c4 c4 c4 c4 c4];
c5=cos(t5);
c5=[c5 c5];
c6=cos(t6);

% Random frequency hops to form a spread signal
spread_signal=[];
spread=[];
for n=1:N
    c=randint(1,1,[1 6]);
    spread=[spread c];
    switch(c)
        case(1)
            spread_signal=[spread_signal c1];
        case(2)
            spread_signal=[spread_signal c2];
        case(3)
            spread_signal=[spread_signal c3];
        case(4)
            spread_signal=[spread_signal c4];
        case(5)        
            spread_signal=[spread_signal c5];
        case(6)
            spread_signal=[spread_signal c6];
    end
end
display('Spreading Code :');
display(spread);
subplot(5,1,4);
plot(1:N*T,spread_signal);
set(gca,'XTick',0:T:N*T,'XLim',[-60 N*T+60],'YLim',[-1.5 1.5]);
title('Spreading Code with 6 frequencies');
ylabel('Amplitude(V)');
xlabel('Time(s)');

% Spreading BPSK Signal into wider band with total of 12 frequencies
freq_hopped_sig=bpsk_sig.*spread_signal;
subplot(5,1,5);
plot(1:N*T,freq_hopped_sig);
set(gca,'XTick',0:T:N*T,'XLim',[-60 N*T+60],'YLim',[-1.5 1.5]);
title('Frequency Hopped Spread Spectrum Signal');
ylabel('Amplitude(V)');
xlabel('Time(s)');

% Receiver Side
figure;
subplot(3,1,1);
plot(1:N*T,freq_hopped_sig);
set(gca,'XTick',0:T:N*T,'XLim',[-60 N*T+60],'YLim',[-1.5 1.5]);
title('Frequency Hopped Spread Spectrum Signal');
ylabel('Amplitude(V)');
xlabel('Time(s)');

% Demodulated BPSK Signal
demod_psk=freq_hopped_sig.*spread_signal;
subplot(3,1,2);
plot(demod_psk);
set(gca,'XTick',0:T:N*T,'XLim',[-60 N*T+60],'YLim',[-1.5 1.5]);
title('Demodulated BPSK Signal');
ylabel('Amplitude(V)');
xlabel('Time(s)');

% Demodulated Binary Signal
demod_sig=[];
for j=0:T:N*T-1
     if demod_psk(j+1)<0
         sig=-ones(1,T);
     else
         sig=ones(1,T);
     end
     demod_sig=[demod_sig sig]; 
end
subplot(3,1,3);
stairs(demod_sig);
set(gca,'XTick',0:T:N*T,'XLim',[-60 N*T+60],'YLim',[-1.5 1.5]);      %Setting axis limits and scale for the graph
title('Demodulated Binary Signal');
ylabel('Amplitude(V)');
xlabel('Time(s)');