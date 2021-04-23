function [fidout1,fidout2] = remove_after_pulse(fileName1,fileName2,dt,bin_size)
%This function takes in two filenames, opens and reads them
%and removes afterpulses that occur at some time difference dt
t_offset=1.5e6 - bin_size;
%read files
%dt=1e5;
fid1=fopen(fileName1,'r');
fid2=fopen(fileName2,'r');
allchannel1=fread(fid1,inf,'uint64');
allchannel2=fread(fid2,inf,'uint64');%-t_offset;

%afterpulse calc
tdif_1=allchannel1(2:length(allchannel1)) - allchannel1(1:(length(allchannel1)-1));
ind1_ap=find(tdif_1>dt);
channel1=allchannel1(ind1_ap);
fileOut1=sprintf('channel1_noafterpulse.dat');
fidout1=fopen(fileOut1,'w+');
fwrite(fidout1,channel1,'uint64');
fseek(fidout1,0,'bof');
tdif_2=allchannel2(2:length(allchannel2)) - allchannel2(1:length(allchannel2)-1);
ind2_ap=find(tdif_2>dt);
channel2=allchannel2(ind2_ap);
fileOut2=sprintf('channel2_noafterpulse.dat');
fidout2=fopen(fileOut2,'w+');
fwrite(fidout2,channel2,'uint64');
fseek(fidout2,0,'bof');
end

