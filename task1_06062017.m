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
t=(0:Tsamp:Tsym-Tsamp); 
alpha=3/100; 
tot_data_bits=100000;
%% Generated data
data=round(rand(1,tot_data_bits));

basic_pulse = sin(2*pi*0.5*fsym*t); % generate basic pulse
basic_pulse= basic_pulse./sqrt(sum(basic_pulse.^2)/fs);
figure (1);
plot(t,basic_pulse);
%% Map to complex signals
%QPSK_bits=-sign(data(1:2:length(data))-0.5)-1i.*sign(data(2:2:length(data))-0.5); %modulated bits
%QPSK_bits= sign(data(1:2:length(data))-0.5)+1i.*sign(data(2:2:length(data))-0.5); %modulated bits
data_one_minus_one=data*2-1; % Map 0 to -1, 1 to 1
I_data=data_one_minus_one(1:2:end-1);
Q_data=data_one_minus_one(2:2:end);
QPSK_bits= -I_data-1i*Q_data; % Map to symbol

%% Build QPSK_pulse_train
QPSK_pulse=basic_pulse'*QPSK_bits;
QPSK_pulse_train=reshape(QPSK_pulse, 1, []);
% figure (2); plot(real(QPSK_pulse_train));
%% Carrier modulation
Modbitslength=length(QPSK_pulse_train);
Msglength=length(data)/2;
Realpart=real(QPSK_pulse_train); 
Imagpart=imag(QPSK_pulse_train);
t1=(0:Tsamp:101*Msglength*Tsamp-Tsamp); 
%% Up conversion
% c=cos(2*pi*fc*t1);
% s=-sin(2*pi*fc*t1);
I=Realpart.*cos(2*pi*fc*t1); 
Q=Imagpart.*-sin(2*pi*fc*t1);
Y=sqrt(2).*(I+Q); %transmitted complex signal
%% Channel build
%r=1/(sqrt(1+alpha.^2));
SNR_dB=0:1:10;    
SNR_lin=10.^(SNR_dB/10);
No = 1./SNR_lin;
No=[-10:1:0];
Ch = zeros(1,length(No));
BER= zeros(1,length(No));
Ch(1) = 1/(sqrt(1+alpha.^2));
Ch(end) = alpha*(1/(sqrt(1+alpha.^2)));
%Ch=[(1/(sqrt(1+alpha.^2))) 0 0 0 0 0 0 0 0 0 0 alpha*(1/(sqrt(1+alpha.^2)))];
%Ch=[r 0 0 0 0 0 0 0 0 0 0 alpha*r]
%% ADDING NOISE
j=0;
for CurrNo = No %dB
CurrNo_lin=10.^(CurrNo/10);%
Var=(CurrNo_lin/2)*fs;
Noise=randn(1,Modbitslength).*sqrt(Var);
j=j+1;

convchannel=conv(Ch, Y);
R=convchannel(1:Modbitslength)+Noise; %here we add the noise
%% DEMODULAtion Down conversion

I_rx=R.*cos(2*pi*fc*t1);
Q_rx=R.*(-sin(2*pi*fc*t1));
I_matchedfilter=conv(I_rx, basic_pulse);
Q_matchedfilter=conv(Q_rx, basic_pulse);
figure (2); 
plot(I_matchedfilter);
%sampling
for k=1:tot_data_bits/2
    I_sampled(k)=I_matchedfilter(1,101*k);
    Q_sampled(k)=Q_matchedfilter(1,101*k);
end

%% demapping
Total_sampled=[I_sampled, Q_sampled];
    for i=1:length(I_sampled)
      if I_sampled(i)>0
        v(i)=[0];
      else v(i)=[1];
    end
    if Q_sampled(i)>0
        u(i)=[0];
    else u(i)=[1];
    end
    end
    
%     for k=1:length(I_sampled)   %------------------------------KB WHY ???  did not understand the logic why at 101 why not 51 to take out the middle value of the pulse
%         if (I_sampled(k)>0)
%             DecodedSignal_I(k)= 0;
%         else
%             DecodedSignal_I(k)=1;
%         end
%     end
%     for k=1:length(Q_sampled)   %------------------------------KB WHY ???  did not understand the logic why at 101 why not 51 to take 
%         if (Q_sampled(k)>0)
%             DecodedSignal_Q(k)= 0;
%         else
%             DecodedSignal_Q(k)=1;
%         end
%     end
    
% DecodedSignal = [DecodedSignal_I; DecodedSignal_Q];    
% Decoded_bits = reshape(DecodedSignal,1,numel(DecodedSignal));
estBits=[v;u];
estBits=reshape(estBits,1,[]);

%% detect errors
% 
% [errors,error_rate]=biterr(data,Decoded_bits);
% error(j)=errors;
% error_sum(j)=errors/length(data);

% 
 xorf=xor(data,estBits); %Decoded_bits
 nr_errors=sum(xorf);
BER(j)= nr_errors;
end
%% Plot The BER vs SNR
ber= BER/tot_data_bits;
No=-10:1:0
figure (15)
semilogy(-No,ber);
title('BER vs SNR') 
ylabel('BER') % x-axis label
xlabel('Eb/No [dB]') % y-axis label

