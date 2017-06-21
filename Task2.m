% Task 2 

clear all;
close all;

%% System variables 
load Signal3.mat
fs= 44100;
fc= 4000;
Tsym= 2.2676e-3;   
fsym=1/Tsym; 
Tsamp=1/fs;
t=(0:Tsamp:Tsym-Tsamp); 
tot_data_bits=100000;
%% Base Pulse 

base_pulse = sin(2*pi*0.5*fsym*t); 
Es=sum(abs(base_pulse).^2)*(1/fs);
P_norm=base_pulse/sqrt(Es);
figure (1);
plot(t,P_norm);

%% Down conversion

t= 0:Tsamp:(length(R)-1)*Tsamp;
Down_I=R.*cos(2*pi*fc*t);
Down_Q=R.*(-sin(2*pi*fc*t));

%% Matched filter 

Match_I=conv(Down_I,P_norm);
Match_Q=conv(Down_Q,P_norm);
plot(sqrt(Match_I.^2+Match_Q.^2));

Signal = Match_I+i*Match_Q;
[Peak start_pos]=max(abs(Signal));

I_sampled=Match_I(start_pos:100:139255);
Q_sampled=Match_Q(start_pos:100:139255);

Signal_sampled= I_sampled+i*Q_sampled;
alpha= Signal_sampled(1)/(2+2*1i);
Signal_Transmitted= Signal_sampled(2:end-1)/alpha;

I_part= real(Signal_Transmitted);
Q_part= imag(Signal_Transmitted);


for k=1:length(I_part)   
    if (I_part(k)>0)
        DemodSignal_I(k)= 0;
    else
        DemodSignal_I(k)=1;
    end
end
for k=1:length(Q_part)  
    if (Q_part(k)>0)
        DemodSignal_Q(k)= 0;
    else
        DemodSignal_Q(k)=1;
    end
end
    
DecmodSignal = [DemodSignal_I; DemodSignal_Q];    
Demod_bits = reshape(DecmodSignal,1,numel(DecmodSignal));

%ASCII%
wh=2.^[6:-1:0];
m=char(Demod_bits(1:7)*wh');
for l=2:floor(length(Demod_bits)/7),
m=[m char(Demod_bits(7*(l-1)+1:7*l)*wh')];
end 
m