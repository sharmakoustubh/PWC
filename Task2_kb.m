% Task 1 
clc;
clear all;
close all;

%% System variables 
load Signal5.mat
fs= 44100;
fc= 4000;
Tsym= 2.2676e-3;  %Symbol time 
fsym=1/Tsym; 
Tsamp=1/fs;
u=0.00023;
t=(0:Tsamp:Tsym-Tsamp); 
alpha=6/100; 
tot_data_bits=100000;

%% Base Pulse yahan se naam change karne hai 

data=round(rand(1,tot_data_bits));
base_pulse = sin(2*pi*0.5*fsym*t); % generate basic pulse
Es=sum(abs(base_pulse).^2)*(1/fs);
Pnorm=base_pulse/sqrt(Es);
figure (1);
plot(t,Pnorm);

%% Down conversion

t= 0:Tsamp:(length(R)-1)*Tsamp;
Down_I=R.*cos(2*pi*fc*t);
Down_Q=R.*(-sin(2*pi*fc*t));

%% matched filter implementation

Match_I=conv(Down_I,Pnorm);
Match_Q=conv(Down_Q,Pnorm);
plot(sqrt(Match_I.^2+Match_Q.^2)); %plot the envelope

Signal = Match_I+i*Match_Q;
[Peak start_pos]=max(abs(Signal));

I_sampled=Match_I(start_pos:100:52920);
Q_sampled=Match_Q(start_pos:100:52920);

Signal_sampled= I_sampled+i*Q_sampled;
alpha= Signal_sampled(1)/(2+2*1i);
Signal_Transmitted= Signal_sampled(2:end-1)/alpha;

I_part= real(Signal_Transmitted);
Q_part= imag(Signal_Transmitted);


for k=1:length(I_part)   
    if (I_part(k)>0)
        DecodedSignal_I(k)= 0;
    else
        DecodedSignal_I(k)=1;
    end
end
for k=1:length(Q_part)  
    if (Q_part(k)>0)
        DecodedSignal_Q(k)= 0;
    else
        DecodedSignal_Q(k)=1;
    end
end
    
DecodedSignal = [DecodedSignal_I; DecodedSignal_Q];    
Decoded_bits = reshape(DecodedSignal,1,numel(DecodedSignal));

%ASCII%
wh=2.^[6:-1:0];
m=char(Decoded_bits(1:7)*wh');
for l=2:floor(length(Decoded_bits)/7),
m=[m char(Decoded_bits(7*(l-1)+1:7*l)*wh')];
end 

m
