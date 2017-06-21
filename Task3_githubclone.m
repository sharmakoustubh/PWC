%clc
close all
clear all
load('signal3.mat')

%% System parameters

fs = 44100;
fc = 10000;
Tsym = 58e-3;
Nsc=128; 
Ncp=20; 
%% Downconversion of the signal

Down_I= sqrt(2).*(R.*cos(2*pi*fc*t));
Down_Q=sqrt(2).*(-R.*sin(2*pi*fc*t));

%% Digital to Analog conversions

[B, A]=butter(7,0.05); 
FilteredSignal_I=filter(B, A, Down_I);
FilteredSignal_Q=filter(B, A, Down_Q);
FilteredSignal =(FilteredSignal_I  +(1i.*FilteredSignal_Q));
figure()
plot(abs(FilteredSignal));

%% Sampling the data

Ts = round((Tsym/Nsc)*fs);
SampledSignal =FilteredSignal(1:Ts:end); 

for t1=1:length(SampledSignal)-(Nsc); 
    m=SampledSignal(t1:t1+(Nsc/2));
    n=SampledSignal(t1+(Nsc/2):t1+Nsc);
    n1=conj(n);
    delay(t1)=abs(sum(m.*n1)/(sqrt(sum((abs(m).^2)))*sqrt(sum((abs(n1).^2)))));
end

%% Removing cyclic prefix and FFT

[Corr start_pos]= max(delay);
PilotSignal=SampledSignal(start_pos:start_pos+Nsc-1);
OFDM_blocks = round(length(SampledSignal(start_pos+Nsc:end))/(Nsc+Ncp))-1;

for i=1: OFDM_blocks 
   OFDM_chunk = SampledSignal(start_pos+(Nsc+Ncp)*i:start_pos+(Nsc+Ncp)*i+Nsc-1);
   FFT_Signal(i,:) = fft(OFDM_chunk);
end

FFT_Pilot=fft(PilotSignal); 
%% generating known pilot

x=zeros(1,Nsc); 
randn('state',100);
P=sign(randn(1,Nsc/2));
x2=2*P; 
x(1:2:end)=x2;
Channel1=FFT_Pilot(1:2:end)./x2;
Channel2=interp1((1:2:Nsc),Channel1,1:Nsc);
%% remove the effect of the channel

for i=1:size(FFT_Signal,1)
     QPSK_Symbols(i,:) = FFT_Signal(i,:)./Channel2; 
end


Signal_Transmitted=(reshape(QPSK_Symbols',[],1))'; 
I_part=sign(real(Signal_Transmitted)); 
Q_part=sign(imag(Signal_Transmitted));

for k=1:length(I_part)   
    if (I_part(k)>0)
        m(k)= 0;
    else
        m(k)=1;
    end
end
for k=1:length(Q_part)  
    if (Q_part(k)>0)
        n(k)= 0;
    else
        n(k)=1;
    end
end
    
Demod_bits = [m; n];    
Demod_bits_stream = reshape(Demod_bits,1,numel(Demod_bits));
%% Viterbi decoding

const_length=6; 
Trellis =poly2trellis(const_length,[77,45]); 
TB_length=6*5;
decoded_bits=vitdec(Demod_bits_stream,Trellis,TB_length,'term','hard'); 

%% Decoding the message

Msg_length=bi2de(decoded_bits(1:10));
Msg_bits=decoded_bits(Nsc+1:end); 
wh=2.^[6:-1:0];
m=char(Msg_bits(1:7)*wh');
for l=2:Msg_length; 
m=[m char(Msg_bits(7*(l-1)+1:7*l)*wh')];
end 
m
