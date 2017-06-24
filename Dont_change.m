clc;
close all;
clear all;
authors='Rajeswari & Jeena & Ardiana';
%% Variables
Fs=44100; % sampling frequency
fc= 5000; % carrier frequency
Nframes= 10; %Number of frames to be transmitted
N_sc=128;  %Number of subcarriers
N_cp=20;  % Number of cyclic prefix
Tsymbol=0.058; % 1 OFDM symbol period
T_sc=Tsymbol/N_sc; % period of 1 subcarrier
F_sc=1/T_sc; % carrier spacing 
deltaf=1/Tsymbol;Rs_ofdm=deltaf*N_sc; %sampling rate or Bandwidth
Nframes=10;
bit_len= Nframes*N_sc;
N=N_sc+N_cp;
%% Generation of bits

% IMAGE TO BITS CONVERSION
I2=imread('NW.jpg'); 
Image = im2uint8(I2,'indexed');
transmit_bits=reshape(de2bi(Image(:)),1,[]);
imshow(Image)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BIT GENERATION
% len = 1000;
% transmit_bits=randi([0 1],1,bit_len); % generating random bits
% coded_bits=transmit_bits;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONVOLUTIONAL ENCODER
trellis=poly2trellis(6,[77,45]); % 6= 5 delays + 1 i/p, 77&45 means 2 o/p
coded_bits=convenc(transmit_bits,trellis);% 77=111 111 means all 5 bits affect 1st o/p, 45=100 101 means 2 delay dont affect 2 o/p
%% Zero Padding
length_zero=rem(length(coded_bits),N_sc*2);%length of zeros to make nultiple of N_sc
num_zero=2*N_sc-length_zero;%number of zeros to be added
zero_padd=zeros(1,round(num_zero));%% generation of zeros
coded_bits=[coded_bits zero_padd];% zero padding/appending
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QPSK MODULATION
qpsk_bits=zeros(1,length(coded_bits)/2);
for k=1:length(qpsk_bits)
    if(coded_bits(2*k-1)==0&&coded_bits(2*k)==0)
        qpsk_bits(k)=1+1i;   
    end
    if(coded_bits(2*k-1)==0&&coded_bits(2*k)==1)
        qpsk_bits(k)=1-1i;  
    end
    if(coded_bits(2*k-1)==1&&coded_bits(2*k)==0)
        qpsk_bits(k)=-1+1i ;  
    end
    if(coded_bits(2*k-1)==1&&coded_bits(2*k)==1)
        qpsk_bits(k)=-1-1i;
    end         
end 
QPSK=qpsk_bits;
%scatterplot(QPSK);grid on;title('qpsk modulated transmitted data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SERIAL TO PARALLEL CONVERSION
% Partitioning of message to blocks of N_sc symbols
% Each subcarrie is modulated as if it was an individual channel
size_message=length(QPSK); % 224128
block_size=floor(size_message/N_sc) %  blocks of msg/channels --> 1751
% Array of zeros with 128 rows and 1751 columns
QPSK_ifft=zeros(N_sc,block_size);
for n=1:block_size % serial to parallel and IFFT
    Ser2Par_signal=QPSK((1:N_sc)+ N_sc *(n-1)); 
    % modulate each subchannel onto the appropiate carrier
    QPSK_ifft(:,n)=ifft(Ser2Par_signal.');  
end
%% CP ADDITION
% A cyclic prefix a repetition of the first section of a symbol that is
% appended to the end of the symbol
OFDMmsg=zeros(N_sc+N_cp,block_size); % 148x1751
for hh=1:block_size % loop through the block size "individual channels"
    OFDMmsg(:,hh)=[QPSK_ifft(end-N_cp+1:end,hh); QPSK_ifft(:,hh)]; % CP+Sym
end
% figure(3)
% plot(real(OFDMmsg))
% title('OFDMmsg with CP Serial/Parallel');grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARALLEL TO SERIAL CONVERSION
% Suming all sub carriers and combining them into one signal
Par2ser_signal=reshape(OFDMmsg,1,[]);%% serial message --> 148x1751
%figure(4);plot(real(Par2ser_signal));title('OFDMmsg with CP');grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADDITION OF PILOTS
% GENERATION OF PILOTS
randn('state',100);
P=sign(randn(1,N_sc/2)); %1x64
x=zeros(1,N_sc);
u=2*P;
x(1:2:end)=2*P; %Every second symbol equals 0 (1x128 pilot symbols)
pilot_ifft=ifft(x); % 1x128
CP_pilot=pilot_ifft(end-N_cp+1:end);
pilot_CP=[CP_pilot pilot_ifft];
%figure(5);plot(real(pilot_CP));title('PilotSignal');grid on;
msg_plus_pilot=[ pilot_CP  Par2ser_signal]; % appending qpsk bits after the pilots
%figure(6);plot(abs(msg_plus_pilot));title('OFDMmsg plus PilotSignal');grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DIGITAL TO ANALOG CONVERSION
%samples/sec*symbols/sec= samples/OFDM symbol
Nsamples= round(Fs/Rs_ofdm); % greatest integer
Interpolated_signal= interp(msg_plus_pilot,Nsamples);
[A,B]=butter(8,(Rs_ofdm/Fs));
%  DAC_signal=Interpolated_signal;
 DAC_signal=filter(A,B,Interpolated_signal);
%figure(7);plot(abs(DAC_signal));title('DAC Signal');grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UP CONVERSION
% t=0:T_sc:T_sc*length(DAC_signal)-T_sc;
t=(0:length(DAC_signal)-1)/Fs;
I=cos(2*pi*fc*t).*sqrt(2);   %I carrier
Q=sin(2*pi*fc*t).*sqrt(2);   %Q carrier
BI=real(DAC_signal).*I;
BQ=imag(DAC_signal).*(-Q);
%OFDM_signal=BI+1i*BQ;
OFDM_signal=BI+1i*BQ;
% figure(8);plot(t,real(OFDM_signal)); xlabel('Time'); ylabel('Amplitude');title('OFDM Signal');grid on;
% figure(9);psd(real(OFDM_signal));title('POWER SPECTRUM');grid on;
R=real(OFDM_signal);
%% SPEAKER
%wavplay(OFDM_signal,Fs);
p = audioplayer(OFDM_signal,fc);
play(p, [1 (get(p, 'SampleRate') * 3)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RECEIVER
% rx_bits=OFDM_rxr(R,t);
% length=length(coded_bits)-length(rx_bits)
% nerr=nnz(coded_bits(:)-rx_bits(:))
% bits2image=bi2de(rx_bits);
% R_I = mat2gray(bits2image);
% imshow(R_I);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOWN CONVERSION
II=cos(2*pi*fc*t).*sqrt(2);   %I carrier
QQ=sin(2*pi*fc*t).*sqrt(2);   %Q carrier
rI=real(OFDM_signal).*II;
rQ=imag(OFDM_signal).*(-QQ);
rr=rI+1i*rQ;
%figure(10);plot(t,rr,'r');grid on;title(['piwc04.m by ',authors])
%% Low pass Butterworth filter
[A1,B1]=butter(9,(Rs_ofdm/Fs),'low'); % Signal_processing_toolbox 
Filtered_r=filter(A1,B1,rr);
% y=Filtered_r; % keeping original signal z
%figure(11);plot(t,Filtered_r,'r');grid on; title(['piwc04.m by ',authors])
%% Sampling
%samples/sec*symbols/sec= samples/OFDM symbol
Nsamples= round(Fs/Rs_ofdm); % greatest integer 
Y=Filtered_r(1:Nsamples:end); %sampling with spacing Tc
t1=t(1:Nsamples:end);
%% Synchronisation
% To find T0
T_periodic= N_sc/2;%the transmitted signal is periodic with period N/2.
% instead of taw in formula put taw=T_periodic-1 and t= T0
% Computing Mye(t) and Normalizing it to 1
myu=zeros(1,(length(Y)-(2*T_periodic)));
for T0=1:(length(Y)-2*T_periodic)
myu(T0)=sum(Y(T0:T0+T_periodic-1)*(Y(T0+T_periodic:T0+2*T_periodic-1))')/(sqrt(sum(abs(Y(T0:T0+T_periodic-1).^2)))*sqrt(sum(abs(Y(T0+T_periodic:T0+2*T_periodic-1).^2))));
end
% close all
%figure(12);plot(abs(myu));
pilot1=find(myu==max(myu))+1  %Pilot1

%% Removal of  Cyclic Prefix
N_OFDM=floor(length(Y(pilot1+N_sc:end))/N);%finding number of OFDM blocks(signal+CP) without pilot,( first prefix from signal)
pilot=Y(pilot1:pilot1+N_sc-1); % pilot_symbol(1:128) from sampled signal
Sig_OFDM=zeros(N_OFDM,N_sc);
Sig_OFDM_fft=zeros(N_OFDM,N_sc);
for n=1:N_OFDM %take away pilot,CP && FFT
    Sig_OFDM(n,:)=Y(pilot1+N*n:pilot1+N*n+N_sc-1);% avoiding pilot+cp : every 128 symbols after every(128 pilot+ 20cp)
   %FFT OF SGNAL
    Sig_OFDM_fft(n,:)=fft(Sig_OFDM(n,:)); %32x128
end
pilot_FFT=fft(pilot);   %1x128
Sig_cp=Sig_OFDM((N_cp+1):end);% 20 first positions of CP shall be removed (signal without CP)
% figure(13);plot(abs(Sig_cp));%without cp
% title('SIgnal without CP');grid on;

 %% Pilot & Channel estimation
%generating the known pilot
randn('state',100);
P1=sign(randn(1,N_sc/2)); %1x64
x1=zeros(1,N_sc);
x1(1:2:end)=2*P1; %Every second symbol equals 0(1x128 pilot symbols)
x2=2*P1;
%%%%%%%%%%%%%%%%%%%%
%% Generating alpha factor

alpha=pilot_FFT(1:2:end)./x2; %% EQUALIZATION
alpha_interpol=interp1((1:2:N_sc),alpha,1:N_sc,'linear','extrap');%1-D interpolation(estimated based on interpolation)
alpha_interpol(128)=alpha(64); %128 cant be estimated from interpolation because it samples  only to 127
% figure(14);
% plot(1:2:length(x1),abs(alpha),'*m','linewidth',2); % Every second symbol will have alpha factors 
% hold on;
% plot(1:length(x1),abs(alpha_interpol),'-g','linewidth',2);
% grid on;
% title(['piwc04.m by ',authors])
%%%%%%%%%%%%%%%%%%%%%
%% DEMODULATION
QPSK_Receiver=zeros(N_sc,N_OFDM);
for n=1:N_OFDM 
QPSK_Receiver(:,n)=Sig_OFDM_fft(n,:)./alpha_interpol;
end
y_re=real(QPSK_Receiver);
y_im=imag(QPSK_Receiver);
MatrixSymb(find(y_re < 0 & y_im < 0))  = 3+1;
MatrixSymb(find(y_re >= 0 & y_im > 0)) = 0+1;
MatrixSymb(find(y_re < 0 & y_im >= 0)) = 2+1;
MatrixSymb(find(y_re >= 0 & y_im < 0)) = 1+1;
% y_hat = zeros(1,length(MatrixSymb));
Matrix=MatrixSymb;
map=[0 0; 0 1; 1 0; 1 1];%% Mapping bits
y2d=map(Matrix,:);
y_hat = reshape(y2d',[],1);
%  rx_data=y_hat;
%% removing zero padding
length_zeroq=rem(length(y_hat),N_sc*2);%length of zeros to make nultiple of N_sc
num_zeros=2*N_sc-length_zeroq;%number of zeros to be added
rx_data=y_hat(1:end-num_zeros);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trellis code(Viterbi Algorithm)(vitdec.m is used)
trelliscode=poly2trellis(6,[77 45]); % constraint length=6 =( 5 delay states+1) .(77,45) rate 1/2 convolutional encoder
tblen=6*5;
%Assumed encoder starts and ends in the all-zero state
rx_bits_int=vitdec(rx_data',trelliscode,tblen,'term','hard');%terminate at zero state-> term,binary input(1's & 0's) ->hard
% rx_bits= uint8(dec2bin(rx_bits_int)');
dff_tx_rx=size(transmit_bits,2)-size(rx_bits_int,2)
nerr=size(find(transmit_bits~=rx_bits_int),2)
%%%%%%%%% 
% %% Image display
% ll = reshape(rx_bits_int, size(rx_bits_int,2)/8, 8);
% xx = bi2de(ll);
% im=reshape(xx,168,150,3);
% imshow(im)