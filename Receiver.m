function RecData = Receiver(Data)
% Ardiana,Jeena & Raji

%% Variables

Fs=44100; % sampling frequency
fc= 5000; % carrier frequency
Nframes= 10; %Number of frames to be transmitted
N_sc=128;  %Number of subcarriers
N_cp=20;  % Number of cyclic prefix
Tsymbol=0.058; % 1 OFDM symbol period
T_sc=Tsymbol/N_sc; % period of 1 subcarrier
deltaf=1/Tsymbol;Rs_ofdm=deltaf*N_sc; %sampling rate or Bandwidth
deltaf=1/Tsymbol;
Rs_ofdm=deltaf*N_sc; %sampling rate or Bandwidth
N=N_sc+N_cp;
time = 10;

%% Recording a message
r=wavrecord(time*Fs,Fs,1); % record time*Fs samples of an audio signal, sampled at rate of Fs. 1 is for default(Stereo,mono)
r=r';
t=(1:length(r))/Fs; % The time of the record

%% DOWN CONVERSION
II=cos(2*pi*fc*t).*sqrt(2);   %I carrier
QQ=sin(2*pi*fc*t).*sqrt(2);   %Q carrier
rI=real(r).*II;
rQ=imag(r).*(-QQ);
rr=rI+1i*rQ;

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

%% trellis code(Viterbi Algorithm)(vitdec.m is used)
trelliscode=poly2trellis(6,[77 45]); % constraint length=6 =( 5 delay states+1) .(77,45) rate 1/2 convolutional encoder
tblen=6*5;
%Assumed encoder starts and ends in the all-zero state
rx_bits_int=vitdec(rx_data',trelliscode,tblen,'term','hard');%terminate at zero state-> term,binary input(1's & 0's) ->hard
rx_bits_int_final = rx_bits_int(1:Data+16);

% rx_bits= uint8(dec2bin(rx_bits_int)');
% dff_tx_rx=size(transmit_bits,2)-size(rx_bits_int,2)
% nerr=size(find(transmit_bits~=rx_bits_int),2)

%% CRC decoding
% creates a CRC detector object, H. This object has the Polynomial property set to POLY, and the other specified properties set to the specified values
% Specify the number of checksums available at each input frame
H = comm.CRCDetector([16 15 2 0],'CheckSumsPerFrame',1);  
[RecData1,CRC]=step(H,rx_bits_int_final');
RecData=RecData1';   

end