%% time stamping code
% Clark Nov to Jan 2020

clear all;
tic

%hist settings
t_max=1e8; %time in ps=1ms for searching function
bin_size=5000; %in ps
t_offset=1.5e6 - bin_size;
bins=round(t_max/bin_size);
histrange=linspace(-t_max,t_max,bins);
co=zeros(1,length(histrange)-1);
histedge=250e3; %time in 500e3ps=500ns for hist
binedge=50;
bins_edge=round(histedge/binedge);
edges=linspace(-histedge,histedge,bins_edge);
N=zeros(1,length(edges)-1);


%get file parameters such as count rate, number of events
%% read in files
fileName1='Results_TimestampsE5_3_4mW_MaxPol_C2.bin';
fileName2='Results_TimestampsE5_3_4mW_MaxPol_C3.bin';
%channel 1
fid1 = fopen(fileName1,'r+');
fseek(fid1,-8,'eof');%go to end of file
last1=fread(fid1,8,'uint64');
events1= ftell(fid1)/8; %total events
%time1=last1/1e12; %total time
rate1=events1/last1; %channel rate
t1=(1/rate1); %average spacing of events in channel
tlast1=last1-t1; %the final event in channel
%channel 2
fid2=fopen(fileName2,'r+');
fseek(fid2,-8,'eof');
last2=fread(fid2,8,'uint64');
events2=ftell(fid2)/8; %total events
%time2=last2/1e12; %total time
rate2=events2/last2; %channel rate
t2=(1/rate2); %average spacing of events in channel
tlast2=last2-t2; %the final event in channel
tdif=tlast2-tlast1; %the difference in time between final events in each channel

%N_ext=rate1*rate2*(bin_size)*(last1);
%N_events=N_ext*bins;

%initialise arrays
%specify number of channels
%N=2;
%specify chunk length to read in to workspace from file
chunk_factor=5e4; %vary this instead of partial.m
chunk1=(events1/chunk_factor);
chunk2=(events2/chunk_factor);
precision = 'uint64';
chunk_plot=chunk_factor*1e2;
dt_size = 0;

%return to beginning of file for reading chunks
fseek(fid1,0,'bof');
fseek(fid2,0,'bof');

for i1=1:chunk_factor
    %read in chunks of each file
    d_array1=fread(fid1,chunk1,'uint64');
    d_array2=fread(fid2,chunk2,'uint64') - t_offset;
    
    %do correlation calculation
    %% find correlations
    
    %cut down channel 1 by removing first and last t_max events
    st=min(d_array1) + t_max;
    fin=max(d_array1) - t_max;
    indices=find(d_array1>st & d_array1<fin);
    ind=length(indices);
    %initialise arrays
    dt=zeros(ind,50); %preallocate dt for speed
    %WHY DOESNT NANS WORK?
    %dtsize=zeros(ind,1); %this will record length of dt i.e number of correlations
    %mini=zeros(ind,1); %this will store the min dt in every loop
    
    for i2=1:ind
        here=indices(i2); %find index value of this iteration
        now=d_array1(here); %find time value of this index
        hi=now + (2*t_max); %set upper bound
        lo=now - (2*t_max); %set lower bound
        %search in next channel using find
        match = find(d_array2<hi & d_array2>lo);
        difference=d_array2(match) - now;
        if isempty(difference) == false
            dt(i2,1:length(difference))=difference;
            %dtsize(i2)=length(difference); %size of dt array
            %dt_size = dt_size + dtsize(i2); %continually add to counter
            %mini(i2)=min(abs(difference)); %minimum dt
        end
    end
    
    %remove nans
    dt_plt=dt(dt~=0);
    %dt_plt = dt(~isnan(dt));
    
    %count data
    [n, edges]= histcounts(dt_plt,edges);%create histogram of pairs file- this command is a bottleneck on program speed, non linear with 'edges' size
    N= N+n;
    [counts,histrange]=histcounts(dt_plt,histrange);
    co=co+counts;
    
    %write data to files
    %     fileDir= 'C:\Users\rache\Documents\Uni\MATLAB\files';
    %     fullfileName=fullfile(fileDir,'dt.dat');
    %     dlmwrite(fullfileName,dt,'-append');
    %     fullfileName=fullfile(fileDir,'dt_size.dat');
    %     dlmwrite(fullfileName,dtsize,'-append');
    %     fullfileName=fullfile(fileDir,'dt_min.dat');
    %     dlmwrite(fullfileName,mini,'-append');
    
    %run('histplot.m');
    %update plot
    i3=mod(i1,round(chunk_plot/chunk_factor));
    if i3==0%only update plot every time chunk_plot_size has been processed- faster
        run('histplot.m');
        
    end
    
    grad2 = ftell(fid2)/max(d_array2);%gradient of fid2_file_pos at end of this chunk/channel2(time)
    new_posn = 8*round(max(d_array1)*grad2/8);%each timestamp has 8 bytes so find time to nearest 8 units
    fseek(fid2,new_posn,'bof');%move file current position to new_posn_B
end


toc
%find the chunk of channel 2 that starts at the same time as channel 1
