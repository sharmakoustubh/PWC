% Task 1 
clc;
clear all;
close all;

%% System variables 

fs=44100;
fc=4000;
Tsym=2.3*10^-3; 
fsym=1/Tsym; 
Tsamp=1/fs;
u=0.00023;
t=(0:Tsamp:Tsym-Tsamp); 
alpha=6/100; 
tot_data_bits=100000;

%% Base Pulse

data=round(rand(1,tot_data_bits));
base_pulse = sin(2*pi*0.5*fsym*t); % generate basic pulse
Es=sum(abs(base_pulse).^2)*(1/fs);
Pnorm=base_pulse/sqrt(Es);
figure (1);
plot(t,Pnorm);

%% Map to complex signals

data_one_minus_one=data*2-1; % Map 0 to -1, 1 to 1
I_data=data_one_minus_one(1:2:end-1);
Q_data=data_one_minus_one(2:2:end);
QPSK_bits= -I_data-1i*Q_data; % Map to symbol

%% Modulation

QPSK_pulse=Pnorm'*QPSK_bits;
QPSK_pulse_train=reshape(QPSK_pulse, 1, []);
Modbitslength=length(QPSK_pulse_train);
Msglength=length(data)/2;
Realpart=real(QPSK_pulse_train); 
Imagpart=imag(QPSK_pulse_train);
t1=(0:Tsamp:101*Msglength*Tsamp-Tsamp); 

%% Up conversion

I=Realpart.*cos(2*pi*fc*t1); 
Q=Imagpart.*-sin(2*pi*fc*t1);
Y=sqrt(2).*(I+Q); %transmitted complex signal

%% SNR and channel parameters

SNR_dB=0:1:15;    
SNR_lin=10.^(SNR_dB/10);
No = 1./SNR_lin;
Ch = zeros(1,length(No));
BER= zeros(1,length(No));
Ch(1) = 1/(sqrt(1+alpha.^2));
Ch(end) = alpha*(1/(sqrt(1+alpha.^2)));

%% Transmission and reception

j=0;
for CurrNo = No 

Var=(CurrNo/2)*fs;
Noise=randn(1,Modbitslength).*sqrt(Var);
j=j+1;

convchannel=conv(Ch, Y);
R=convchannel(1:Modbitslength)+Noise; %here we add the noise

%% Demodulation and Down conversion

I_rx=R.*cos(2*pi*fc*t1);
Q_rx=R.*(-sin(2*pi*fc*t1));
I_matchedfilter=conv(I_rx, Pnorm);
Q_matchedfilter=conv(Q_rx, Pnorm);
figure (2); 
plot(I_matchedfilter);

for k=1:tot_data_bits/2
    I_sampled(k)=I_matchedfilter(1,101*k);
    Q_sampled(k)=Q_matchedfilter(1,101*k);
end


Total_sampled=[I_sampled, Q_sampled];

for k=1:length(I_sampled)   
    if (I_sampled(k)>0)
        DecodedSignal_I(k)= 0;
    else
        DecodedSignal_I(k)=1;
    end
end
for k=1:length(Q_sampled)  
    if (Q_sampled(k)>0)
        DecodedSignal_Q(k)= 0;
    else
        DecodedSignal_Q(k)=1;
    end
end
    
DecodedSignal = [DecodedSignal_I; DecodedSignal_Q];    
Decoded_bits = reshape(DecodedSignal,1,numel(DecodedSignal));

%% Error Detection
% 
[errors,error_rate]=biterr(data,Decoded_bits);
error(j)=errors;
error_sum(j)=errors/length(data);

BER(j)= error_rate
end

%% Plot The BER vs SNR
ber= BER/tot_data_bits;
figure (15)
semilogy(SNR_dB,BER); % ber %error_sum
title('BER vs SNR') 
ylabel('BER') % x-axis label
xlabel('Eb/No [dB]') % y-axis label

