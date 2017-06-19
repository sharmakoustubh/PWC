%clc
close all
clear all
load('signal3.mat')

%% System parameters

fs = 44100; % The sapling Frequency
fc = 10000; % The carrier frequency
Tsym = 58e-3; % The symbol Time for 1 OFDM
Nsc=128; % 128 subcarriers
Ncp=20; %  20 the cyclic lenght

% R is the received signal from the signal 3 and t is the time taken from
% signal 3
%% Downconversion of the signal

Down_I= sqrt(2).*(R.*cos(2*pi*fc*t));
Down_Q=sqrt(2).*(-R.*sin(2*pi*fc*t));


%% Digital to Analog conversions
% A low pass filter(Butterworth filter)
%load 'Butter.mat';
[B, A]=butter(7,0.05); % designs a lowpass,Butterworth filter, 7.0 is the filter order - integer scalar, 0.05 is the cutoff frequency and low is the filtertype
%[B, A]=butter(8,0.05);

FilteredSignal_I=filter(B, A, Down_I); % filters the input data, Down_I, using a rational transfer function defined by the numerator and denominator coefficients B and a, respectively.
FilteredSignal_Q=filter(B, A, Down_Q); % filters the input data, Down_I, using a rational transfer function defined by the numerator and denominator coefficients B and a, respectively.
FilteredSignal =(FilteredSignal_I  +(1i.*FilteredSignal_Q));

figure()
plot(abs(FilteredSignal));

%% Sampling the data
% Sampling convert analog to signal, sampling frequency is 20 time larger
% than our frequency of our sample is ok if its before or after. Sampling
% frequency much larger than the data
Ts = ceil(fs/(1/(Tsym/Nsc))); % rounds each element of X to the nearest integer greater than or equal to that element.
sampling_data=FilteredSignal(1:Ts:end); % so there will be 6353

% Computing Mye(t) and Normalizing it to 1
for time=1:length(sampling_data)-(Nsc);
    % To comute mye
    OFDM_sig =(sampling_data(time:time+Nsc/2)*(sampling_data(time+Nsc/2:time+Nsc))');
    mye_sum=sum(OFDM_sig);
    % to normalize to 1
    del = (sqrt(sum(abs(sampling_data(time:time+Nsc/2)).^2))*sqrt(sum(abs(sampling_data(time+Nsc/2:time+2*Nsc/2).^2))));
    delay(time) = mye_sum/del;
end

%% Pilots and removing cyclic prefix

% find the start pilot
start_pilot=find(abs(delay)== max(abs(delay)));% P=6288  find the start of the Pilot

pilot_sig=sampling_data(start_pilot:start_pilot+Nsc-1);% the pilot samples N=128 symbols, if earlier ok
ofdm_first_cypr=sampling_data(start_pilot+Nsc:end); % the first cyclic prefix of the signal

% finding the number of OFDM blocks which is the signal + cyclix lenght
% without the pilot --> first prefix
OFDM_blocks = round(length(ofdm_first_cypr)/(Nsc+Ncp))-1; % Round to the nearest integer

% Removing the pilot the cyclic prefix and the fft
for i=1: OFDM_blocks 
   % to avoid the pilot+cyclic lenght : every 128 symbols after every 128(subcarrier)+ 20(Cycic lenght)
   OFDM_sig = start_pilot+(Nsc+Ncp)*i:start_pilot+(Nsc+Ncp)*i+Nsc-1;
   cp = sampling_data(OFDM_sig);
   fourier_signal(i,:) = fft(cp);
   
end
fourier_pil=fft(pilot_sig); 


%% generating known pilot

% The overhead of the block
% generating the known pilot
x=zeros(1,Nsc); 
randn('state',100);
P=sign(randn(1,Nsc/2));
x2=2*P; % the transmitted signal
x(1:2:end)=x2;
%Half of the sub channel have pilot symbol
% the other half is based on interpolation
channel=fourier_pil(1:2:end)./x2; % Every second pilot is 0 --> 64
channel_2=interp1((1:2:Nsc),channel,1:Nsc); % The Interpolate to have 128

%plot (1:2:128,abs(channel), 'or');

%% The demodulation

for i=1:size(fourier_signal,1)
    % dividing the received message by channel
    de_mod(i,:) = fourier_signal(i,:)./channel_2; 
end

% to create one row matrix
demod=(reshape(de_mod',[],1))'; 

% To convert to -1 and 1
decode_I=sign(real(demod)); 
decode_Q=sign(imag(demod));


% bits received I
bits_rec_I=zeros(1,length(decode_I));
for j=1:length(decode_I)
    if (decode_I(j)==1)
        bits_rec_I(j)=0;
    else
        bits_rec_I(j)=1;
    end
end  

% Bits received Q 
bits_rec_Q=zeros(1,length(decode_Q));
for j=1:length(decode_Q)
    if (decode_Q(j)==1)
        bits_rec_Q(j)=0;
    else
        bits_rec_Q(j)=1;
    end
end 

symbol_bits=[bits_rec_I;bits_rec_Q];

received_bits=reshape(symbol_bits,1,[]); % To get the bits in one row
        
%% Trellis code 

% The Viterbi Algorithm is used (vitdec.m)
% rate 1/2 convolutional encoder
trellscode=poly2trellis(6,[77 45]); % constraint length=6 5 delay states +1
tblen=6*5;
% we assume that the encoder starts and ends in the all-zero state
decoded=vitdec(received_bits,trellscode,tblen,'term','hard'); % terminate at zero state -> binary input(1's & 0's)

%% Decoding the message

length_Message=bi2de(decoded(1:10)); % The lenght of the message
f_bits=decoded(Nsc+1:end); 
wh=2.^[6:-1:0];
Message=char(f_bits(1:7)*wh');
for l=2:length_Message; %floor(length(f_bits)/7),
Message=[Message char(f_bits(7*(l-1)+1:7*l)*wh')];
end 
Message