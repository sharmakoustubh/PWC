clear all;
close all;
warning('off');

tic;
Header = 10;
Data = 25;
Rx_buffer = Header + Data;
m = 0;

Ack = [1 1 1 1 1];
Nack = [1 0 1 0 1];
Picture_data = [];

%% Control Packet Rx & decoding
while(1)
    while(1)
        T_RxSyn = floor(rem(toc, 10));
        if(T_RxSyn == 4)
            disp('Listening for Control Packet...');
            break;
        end
    end
    
   [Rx_bits_CP, Rx_CRC_CP] = Rx_func(15);
      
   if(Rx_CRC_CP == 0)
       disp('Control Packet Correct');
       m = 1;
       str_conf = [num2str(bi2de(Rx_bits_CP(10:end))) ' Packets, ' num2str(bi2de(Rx_bits_CP(1:9))) ' bit per Packet'];
       disp(str_conf);
       pause(1);
       
       while(1)
             T_RxSyn = floor(rem(toc, 10));
             if(T_RxSyn == 9)
                 disp('Sending ACK...');
                 disp(' ');
                 break;
             end
       end
       
       pause(1);
       Tx_func(Ack);
       break;
       
   else
       disp('Control Packet Error');
       pause(2);
       
       while(1)
             T_RxSyn = floor(rem(toc, 10));
             if(T_RxSyn == 9)
                 disp('Sending NACK...');
                 disp(' ');
                 break;
             end
       end
       
       Tx_func(Nack);       
   end   
end

while(1)
    while(1)
        T_RxSyn = floor(rem(toc, 10));
        if(T_RxSyn == 4)
            str = ['Listening for Picture Packet #' num2str(m)];
            disp(str);
            break;
        end
    end
    
    [Rx_bits, Rx_CRC] = Rx_func(Rx_buffer);
    
    %Rx_bits
    
    if(Rx_CRC == 0)
       str2 = ['Packet #' num2str(bi2de(Rx_bits(1:10))) ' correct'];
       disp(str2);
       
       Picture_data = [Picture_data Rx_bits(Header+1:Rx_buffer)];
       
       while(1)
             T_RxSyn = floor(rem(toc, 10));
             if(T_RxSyn == 9)
                 disp('Sending ACK...');
                 disp(' ');
                 break;
             end
       end
        
       Tx_func(Ack);
       m = bi2de(Rx_bits(1:10));
       m = m+1;
       %break;
                 
    else
       str3 = ['Packet #' num2str(m) ' Error'];
       disp(str3);
       pause(2);
       
       while(1)
             T_RxSyn = floor(rem(toc, 10));
             if(T_RxSyn == 8)
                 disp('Sending NACK...');
                 disp(' ');
                 break;
             end
       end
       
       Tx_func(Nack);
    end
end

Picture_Rx_crc = reshape(Picture_data,10,10);
imshow(Picture_Rx_crc)