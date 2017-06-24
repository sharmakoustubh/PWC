function [error,OUTdata,crc_result] = Receiving(packet_length,playtime)
    %load('Butter.mat');
    fs=44100;%Hz sampling Hz
    fc=10e03;%Hz carrier Hz
    Nsc=128;%number of subcarriers in OFDM. Ts is the transmissiontime of Nsc data symbols
    Ncp=20;%length of the cyclic prefix
    N=Nsc+Ncp; %length of OFDM block with cp
    Tofdm=58e-3;%symbol time fo 1 OFDM symbol without cp
    SRofdm=1/(Tofdm/Nsc);%Hz sample rate of the OFDM 2206.9 'bw'
    error = 0;

    %% Record msg
   % R=wavrecord(playtime*fs,fs,1);
    recorder=audiorecorder(fs,16,1);
    recordblocking(recorder,18);
    R=getaudiodata(recorder); 
    %R=wavrecord(8*fs,fs,1);
    R=R';
    t = (1:length(R))/fs;

    %% Down conversion
    cI=cos(2*pi*fc*t).*sqrt(2);   %I carrier
    cQ=sin(2*pi*fc*t).*sqrt(2);   %Q carrier
    rI_Rx=R.*cI;
    rQ_Rx=R.*(-cQ);
    rS_Rx=rI_Rx+1i*rQ_Rx; %1x17760 132300=44100x3

    %plot(t,rS);
    %% LOWPASS FILTER
    [B,A]=butter(8,SRofdm/fs,'low'); %license for Signal_processing_toolbo missing
   
    r_filt=filter(B,A,rS_Rx);

    %% Sample
    samPsym=round(fs/SRofdm);
    r_samp=r_filt(10:samPsym:end);
    
    %% sync
    Tper=Nsc/2;
    miu=zeros(1,length(r_samp)-2*Tper);
    for T=1:(length(r_samp)-2*Tper)
        miu(T)=sum(r_samp(T:T+Tper-1).*conj(r_samp(T+Tper:T+2*Tper-1)))/(sqrt(sum(abs(r_samp(T:T+Tper-1).^2)))*sqrt(sum(abs(r_samp(T+Tper:T+2*Tper-1).^2))));
    end
     
    %% remove cp
    Pilot_position=find(miu==max(miu));
    Nofdm=floor(length(r_samp(Pilot_position-Ncp-1:end))/N);%number of potential OFDM blocks
    if(Nofdm <= 1)
        error = 1;
        return;
    end
    r_ofdm=zeros(Nofdm,Nsc);
    for ii=1:Nofdm
        r_ofdm(ii,:)=r_samp(1,(Pilot_position+N*(ii-1)):(Pilot_position+N*(ii-1)+Nsc-1)); %Nx128
    end
    %% FFT
    r_ofdm_fft=zeros(Nofdm,Nsc);
    for ii=1:Nofdm
        r_ofdm_fft(ii,:)=fft(r_ofdm(ii,:)); %Nx128
    end
    r_ofdm_fft_temp=transpose(r_ofdm_fft); %128xN
    r_msg_fft=reshape(r_ofdm_fft_temp,1,numel(r_ofdm_fft_temp));

    pilot_fft=r_msg_fft(1:Nsc);
    r_msg=r_msg_fft(Nsc+1:end);

    %% Pilot& Channel estimation
    %generate pilot
    randn('state',100);
    P=sign(randn(1,Nsc/2)); %1x64
    x=zeros(1,Nsc);
    x(1:2:end)=2*P; %1x128 pilot symbol
    %estimate channel
    Alpha=pilot_fft(1:2:end)./(2*P);    %1x64./1x64=1x64
    Alpha_interpolation=interp1((1:2:Nsc),Alpha,(1:Nsc)); %1-D interpolation
    Alpha_interpolation(128)=Alpha(64); %??n128 cant be estimate from interpolation cause sample to 127
   
    for i=1:Nofdm-1
        r_qpsk_afterCH(Nsc*(i-1)+1:Nsc*i)=r_msg(Nsc*(i-1)+1:Nsc*i)./Alpha_interpolation; %128x32=4096, 1OFDM symbol=4096 QPSK symbol
    end
    %% QPSK demodulation
    r_qpsk_temp=r_qpsk_afterCH;
    zI=real(r_qpsk_temp);
    zQ=imag(r_qpsk_temp);
    zI(zI>0)=0;
    zI(zI<0)=1;
    zQ(zQ>0)=0;
    zQ(zQ<0)=1;
    r_qpsk_de=[zI;zQ];
    r_VA=reshape(r_qpsk_de,1,numel([zI,zQ]));
    %% VA decoding
    trellis=poly2trellis(6,[77 45]);
    tblen=6*5;
    r_ascii=vitdec(r_VA,trellis,tblen,'term','hard');
    r_ascii_fix=r_ascii(1:packet_length+16);
    %r_ascii_fix=r_ascii(1:7264+16);

    %% CRC decoding
    H2 = comm.CRCDetector([16 15 2 0],'CheckSumsPerFrame',1);
    [r_crc,error_crc]=step(H2,r_ascii_fix');
    r_crc=r_crc';
        if (error_crc==0)
            crc_result=1;
        else
            crc_result=0;
        end
        OUTdata=r_crc;
end