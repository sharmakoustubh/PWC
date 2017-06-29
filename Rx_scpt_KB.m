clear all;
close all;
%warning('off');
tic;
Ack = [1 1 1 1 1];
Nack = [1 0 1 0 1];
Img_bits = [];
Packet_size_info = 10;
data = 25;
Total_packet_size = Packet_size_info + data;
i = 0;

%% Synchronisation
while(1)
    %% receive sync packets-------------
    while(1)
        Time_taken = floor(rem(toc, 10));
        if(Time_taken == 4)
            toc
            fprintf('\nWaiting for Sync Packet \n');
            break;
        end
    end
    [Received_data, Errors] = Rx_func_KB(15);
    if(Errors == 0)
        fprintf('Sync Packet was correct \n');
        i = 1;
        String1 = [num2str(bi2de(Received_data(10:end))) ' Packets \n'];
        fprintf(String1);
        pause(1);
        while(1)
            Time_taken = floor(rem(toc, 10));
            if(Time_taken == 9)
                toc
                fprintf('\nTransmitting ACK \n');
                
                break;
            end
        end
        pause(1);
        Tx_func_KB(Ack);
        break;
    else
        fprintf(' \nSync Packet Error \n');
        pause(2);
        while(1)
            Time_taken = floor(rem(toc, 10));
            if(Time_taken == 9)
                toc
                fprintf('Transmitting NACK \n');
                fprintf(' ');
                break;
            end
        end
        Tx_func_KB(Nack);
    end
end
%% Receiving Image packets
while(1)
    %% receive image data packets-------------
    while(1)
        Time_taken = floor(rem(toc, 10));
        if(Time_taken == 4)
            toc
            String2 = ['Waiting for data Packet: ' num2str(i),'\n'];
            fprintf(String2);
            break;
        end
    end
    [Received_data, Errors] = Rx_func_KB(Total_packet_size);
    %% Transmit ACK or NACK-------------
    if(Errors == 0)
        String3 = ['Packet:' num2str(bi2de(Received_data(1:10))) ' correct \n'];
        fprintf(String3);
        Img_bits = [Img_bits Received_data(Packet_size_info+1:Total_packet_size)];
        while(1)
            Time_taken = floor(rem(toc, 10));
            if(Time_taken == 9)
                toc
                fprintf('Transmitting ACK \n');
                break;
            end
        end
        Tx_func_KB(Ack);
        i = bi2de(Received_data(1:10));
        i = i+1;
        
    else
        String4 = [' \nPacket:' num2str(i) ' Error \n'];
        fprintf(String4);
        pause(2);
        while(1)
            Time_taken = floor(rem(toc, 10));
            if(Time_taken == 8)
                toc
                fprintf('Transmitting NACK \n');
                break;
            end
        end
        Tx_func_KB(Nack);
    end
    if (i==5)
        Image = reshape(Img_bits,10,10);
        imshow(Image)
        break;
    end
    
end
Image = reshape(Img_bits,10,10);
imshow(Image)