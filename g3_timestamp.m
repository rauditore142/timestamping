%% time stamping code
% g3 integration
% Clark Nov to Jan 2020

clear all;
tic

%hist settings
%change these as required- bin size, offset, searching time
t_max=1e6; %time in ps=1ms for searching function
bin_size=5000; %in ps
%t_offset=2.5e6 - bin_size; bunched
t_offset=1.5e6 - bin_size; %antibunched
bins=round(t_max/bin_size);
%full range
histrange=linspace(-t_max,t_max,bins);
N1=zeros(1,length(histrange)-1);
N2=zeros(1,length(histrange)-1);
N3=zeros(length(histrange),length(histrange));
%set up plots
%full range plot
% will have to use surf or contour?
%pick a reference channel - eg C1 - and plot dx=dt(C1-C2), dy=dt(C1-C3)
h=surf(histrange,histrange,N3);


%% read in files
%get file parameters such as count rate, number of events
%enter file name here
fileName1='bunched_test_g4_0.7V_10min_C1.bin';
fileName2='bunched_test_g4_0.7V_10min_C2.bin';
fileName3='bunched_test_g4_0.7V_10min_C3.bin';
%channel 1
fid1 = fopen(fileName1,'r+');
fseek(fid1,-8,'eof');%go to end of file
last1=fread(fid1,8,'uint64');
events1= ftell(fid1)/8; %total events
rate1=events1/last1; %channel rate
t1=(1/rate1); %average spacing of events in channel
tlast1=last1-t1; %the final event in channel
%channel 2
fid2=fopen(fileName2,'r+');
fseek(fid2,-8,'eof');
last2=fread(fid2,8,'uint64');
events2=ftell(fid2)/8; %total events
rate2=events2/last2; %channel rate
t2=(1/rate2); %average spacing of events in channel
tlast2=last2-t2; %the final event in channel
tdif=tlast2-tlast1; %the difference in time between final events in each channel
%channel 3
fid3 = fopen(fileName3,'r+');
fseek(fid3,-8,'eof');%go to end of file
last3=fread(fid3,8,'uint64');
events3= ftell(fid3)/8; %total events
rate3=events3/last3; %channel rate
t3=(1/rate3); %average spacing of events in channel
tlast3=last3-t3; %the final event in channel

%initialise arrays
%specify chunk length to read in to workspace from file
%5e4 recommended
chunk_factor=5e4; 
chunk1=(events1/chunk_factor);
chunk2=(events2/chunk_factor);
chunk3=(events3/chunk_factor);
precision = 'uint64';
chunk_plot=chunk_factor*1e2;
N_exp=zeros(chunk_factor,1);
co_exp=zeros(chunk_factor,1);
t_int=zeros(chunk_factor,1);

%return to beginning of file for reading chunks
fseek(fid1,0,'bof');
fseek(fid2,0,'bof');
fseek(fid3,0,'bof');

for i1=1:chunk_factor
    %read in chunks of each file
    d_array1=fread(fid1,chunk1,'uint64');
    d_array2=fread(fid2,chunk2,'uint64');% - t_offset;
    d_array3=fread(fid3,chunk3,'uint64');
    %% find correlations
    %cut down channel 1 by removing first and last t_max events
    st=min(d_array1) + t_max;
    fin=max(d_array1) - t_max;
    indices=find(d_array1>st & d_array1<fin);
    ind=length(indices);
    %initialise arrays
    dt12_array=zeros(ind,50); %preallocate dt for speed
    dt13_array=zeros(ind,50);
    %WHY DOESNT NANS WORK?
    %dtsize=zeros(ind,1); %this will record length of dt i.e number of correlations
    %mini=zeros(ind,1); %this will store the min dt in every loop
%     t_int(i1)=(d_array1(end)+d_array2(end))/2;
%     N_exp(i1)=rate1*rate2*binedge*t_int(i1);
  
    
    for i2=1:ind
        here=indices(i2); %find index value of this iteration
        now=d_array1(here); %find time value of this index
        hi=now + (2*t_max); %set upper bound
        lo=now - (2*t_max); %set lower bound
        %search in next channel using find
        match2 = find(d_array2<hi & d_array2>lo);
        if isempty(match2) == false
            match3 = find(d_array3<hi & d_array3>lo);
            if isempty(match3)==false
                dt12=d_array2(match2) - now;
                dt13=d_array3(match3) - now;
                dt12_array(i2,1:length(dt12))=dt12;
                dt13_array(i2,1:length(dt13))=dt13;
                %dtsize(i2)=length(dt1); %size of dt array
                %mini(i2)=min(abs(dt1)); %minimum dt
            end 
        end
    end
    
    %remove nans
    dt_plt12=dt12_array(dt12_array~=0);
    dt_plt13=dt13_array(dt13_array~=0);
    %dt_plt = dt(~isnan(dt));
    
    %count data
    [n1, histrange]= histcounts(dt_plt12,histrange);%create histogram of pairs file- this command is a bottleneck on program speed, non linear with 'edges' size
    N1= N1+n1;
    [n2,histrange]=histcounts(dt_plt13,histrange);
    N2=N2+n2;
    
   
    for i3=1:length(histrange)-1 %column iterator
        for i4=1:length(histrange)-1 %row iterator
            pos1=N1(i3);
            pos2=N2(i3);
            if pos1~=0
                if pos2~=0
                    N3(i4,i3)=N1(i3)+N2(i3);
                end
            end
        end
    end
    
    %write data to files
    %     fileDir= 'C:\Users\rache\Documents\Uni\MATLAB\files';
    %     fullfileName=fullfile(fileDir,'dt.dat');
    %     dlmwrite(fullfileName,dt,'-append');
    %     fullfileName=fullfile(fileDir,'dt_size.dat');
    %     dlmwrite(fullfileName,dtsize,'-append');
    %     fullfileName=fullfile(fileDir,'dt_min.dat');
    %     dlmwrite(fullfileName,mini,'-append');
    
    
    %% update plot
    i3=mod(i1,round(chunk_plot/chunk_factor));
    if i3==0%only update plot every time chunk_plot_size has been processed- faster
        h.ZDataSource = 'N3';
        refreshdata;
        shg;
    end
    
    grad2 = ftell(fid2)/max(d_array2);%gradient of fid2_file_pos at end of this chunk/channel2(time)
    new_posn_2 = 8*round(max(d_array1)*grad2/8);%each timestamp has 8 bytes so find time to nearest 8 units
    fseek(fid2,new_posn_2,'bof');%move file current position to new_posn_B
    grad3=ftell(fid3)/max(d_array3);
    new_posn_3=8*round(max(d_array2)*grad2/8);
    fseek(fid3,new_posn_3,'bof');
end


fclose('all');
toc
