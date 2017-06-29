clear all;
close all;
%warning('off');

tic;
%% read image
image = imread('test2.bmp');
Ack = [1 1 1 1 1];
i = 1;
Packet_size_info = 10;
data = 25;
NoOfPackets = ceil(numel(image)/data);
Len_packet = de2bi(data + Packet_size_info,9);
Number_packet = de2bi(NoOfPackets,6);
Sync_packet = [Len_packet Number_packet];

%% Synchronisation
String1 = [' \nTotal Packets to be tranmitted are: ' num2str(NoOfPackets),'\n'];
fprintf(String1);

while(1)
    %% send sync packets-------------
    while(1)
        Time_taken = floor(rem(toc, 10));
        if(Time_taken == 4)
            toc
            fprintf('Transmitting Sync Packet \n');
            break;
        end
    end
    Tx_func_KB(Sync_packet);
    fprintf(' \ndone  \n');
    %% wait for ACK or NACK from Receiver-------------
    while(1)
        Time_taken = floor(rem(toc, 10));
        if(Time_taken == 9)
            toc
            fprintf('\nWaiting for ACK or NACK for Sync Packet \n');
            break;
        end
    end
    [ACKorNACK, Error] = Rx_func_KB(5);
    if(ACKorNACK == Ack)
        fprintf('ACK Received \n');
        break;
    else
        fprintf('\nNACK Received or no response received \n');
    end
end
%% Sending Image packets
while(1)
    %% send data packets-------------
    while(1)
        Data_packet = [de2bi(i,10) image((i-1)*data+1:i*data)];
        Time_taken = floor(rem(toc, 10));
        if(Time_taken == 4)
            toc
            String2 = ['\nTransmitting data Packet:' num2str(i)];
            fprintf(String2);
            break;
        end
    end
    Tx_func_KB(Data_packet);
    String3 = [' \nData Packet ',num2str(i),' sent'];
    fprintf(String3);
    %% wait for ACK or NACK from Receiver-------------
    while(1)
        Time_taken = floor(rem(toc, 10));
        if(Time_taken == 8)
            toc
            String4 = [' \nWaiting for ACK or NACK for data Packet:' num2str(i)];
            fprintf(String4);
            break;
        end
    end
    [Rx_bits, Error] = Rx_func_KB(5);
    if(Rx_bits == Ack)
        fprintf(' \nACK Received');
        fprintf(' ');
        i = i+1;
    else
        fprintf(' \nNACK  or no response received \n');
    end
    if(i == NoOfPackets +1)
        fprintf(' \nImage transferred');
        break;
    end
end

