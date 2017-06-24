clear all;
close all;
warning('off');

tic;
n = 1;

Header = 10;
Data = 25;

%% Image read & conversion
Picture = imread('test2.bmp');
Nr_packets = ceil(numel(Picture)/Data);
Ack = [1 1 1 1 1];

%% Control Packet
L_pack_bin = de2bi(Data + Header,9);
Nr_pack_bin = de2bi(Nr_packets,6);

C_pack = [L_pack_bin Nr_pack_bin];

%% Synchro & trans of Control Packet
str0 = ['Total Packets: ' num2str(Nr_packets) ', Packet length: ' num2str(Data + Header)];
disp(str0);

while(1)
    while(1)
        T_TxSyn = floor(rem(toc, 10));
        if(T_TxSyn == 4)
            disp('Sending Control Packet...');
            break;
        end
    end
    
    Tx_func(C_pack);
    disp('Control Packet sent.');
        
    while(1)
        T_TxSyn = floor(rem(toc, 10));
        if(T_TxSyn == 9)
            disp('Receiving ANS for Control Packet...');
            break;
        end
    end  
  
    [Rx_bits_CP, Rx_CRC_CP] = Rx_func(5);     
    
     if(Rx_bits_CP == Ack)
         disp('ACK Received');
         disp(' ');
         break;
     else
         disp('NACK Received');
         disp(' ');
     end
end

while(1)
    while(1)
        Pic_pack = [de2bi(n,10) Picture((n-1)*Data+1:n*Data)]; 
        T_TxSyn = floor(rem(toc, 10));
        if(T_TxSyn == 4)
            str = ['Sending Picture Packet #' num2str(n)];
            disp(str);
            break;
        end
    end
    
    Tx_func(Pic_pack);
    str2 = ['Picture Packet #',num2str(n),' sent'];
    disp(str2);
    
    while(1)
        T_TxSyn = floor(rem(toc, 10));
        if(T_TxSyn == 8)
            str3 = ['Receiving ANS for Picture Packet #' num2str(n) ' ...'];
            disp(str3);
            break;
        end
    end
    
    [Rx_bits, Rx_CRC] = Rx_func(5);
     
    if(Rx_bits == Ack)
         disp('ACK Received');
         disp(' ');
         n = n+1;               
    else
         disp('NACK Received');
         disp(' ');
    end
    
        
    %if(n == Nr_packets +1)
    %    disp('Picture sent');
    %    break;
    %end
        
    
end

