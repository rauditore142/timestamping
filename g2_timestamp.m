%% time stamping code
% Clark Nov to Jan 2020

clear all;
tic

%hist settings
%change these as required- bin size, offset, searching time
t_max=1e8; %time in ps=1ms for searching function
bin_size=500; %in ps
%t_offset=2.5e6 - bin_size; bunched
t_offset=1.5e6 - bin_size; %antibunched
bins=round(t_max/bin_size);
%full range
histrange=linspace(-t_max,t_max,bins);
co=zeros(1,length(histrange)-1);
%`zoomed range' - comment out if not required
histedge=1e6; %time in 500e3ps=500ns for hist
binedge=50;
bins_edge=round(histedge/binedge);
edges=linspace(-histedge,histedge,bins_edge);
N=zeros(1,length(edges)-1);

%set up plots
%full range plot
s1=subplot(1,2,1);
h1=plot(histrange(1:length(co)),co,'.-');
xlabel('Time in ps');
ylabel('Counts');
%zoomed in plot - comment out if not required
s2=subplot(1,2,2);
h2=plot(edges(1:length(N)),N,'.-');
xlabel('Time in ps');
ylabel('Counts');

%% read in files
%get file parameters such as count rate, number of events
%enter file name here
fileName1='Results_TimestampsE5_3_4mW_MaxPol_C2.bin';
fileName2='Results_TimestampsE5_3_4mW_MaxPol_C3.bin';
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

%for normalisation
N_exp=rate1*rate2*binedge*last1;
co_exp=rate1*rate2*bin_size*last1;

%initialise arrays
%specify chunk length to read in to workspace from file
%5e4 recommended
chunk_factor=5e4; 
chunk1=(events1/chunk_factor);
chunk2=(events2/chunk_factor);
precision = 'uint64';
chunk_plot=chunk_factor*1e2;

%return to beginning of file for reading chunks
fseek(fid1,0,'bof');
fseek(fid2,0,'bof');

for i1=1:chunk_factor
    %read in chunks of each file
    d_array1=fread(fid1,chunk1,'uint64');
    d_array2=fread(fid2,chunk2,'uint64') - t_offset;
    
    %% find correlations
    %cut down channel 1 by removing first and last t_max events
    st=min(d_array1) + t_max;
    fin=max(d_array1) - t_max;
    indices=find(d_array1>st & d_array1<fin);
    ind=length(indices);
    %initialise arrays
    dt=zeros(ind,50); %preallocate dt for speed
    %WHY DOESNT NANS WORK?
    dtsize=zeros(ind,1); %this will record length of dt i.e number of correlations
    mini=zeros(ind,1); %this will store the min dt in every loop
    
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
            dtsize(i2)=length(difference); %size of dt array
            mini(i2)=min(abs(difference)); %minimum dt
        end
    end
    
    %remove nans
    dt_plt=dt(dt~=0);
    %dt_plt = dt(~isnan(dt));
    
    %count data
    [n, edges]= histcounts(dt_plt,edges);%create histogram of pairs file- this command is a bottleneck on program speed, non linear with 'edges' size
    N= N+n;
    N_plot=N/N_exp;
    [counts,histrange]=histcounts(dt_plt,histrange);
    co=co+counts;
    co_plot=co/co_exp;
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
       x1=histrange(1:length(co));
        h1.XDataSource = 'x1';
        h1.YDataSource = 'co_plot';
        refreshdata;
        %update plot 2
        %comment out if not required
        x2=edges(1:length(N));
        h2.XDataSource = 'x2';
        h2.YDataSource= 'N_plot';
        refreshdata;
        shg;
    end
    
    grad2 = ftell(fid2)/max(d_array2);%gradient of fid2_file_pos at end of this chunk/channel2(time)
    new_posn = 8*round(max(d_array1)*grad2/8);%each timestamp has 8 bytes so find time to nearest 8 units
    fseek(fid2,new_posn,'bof');%move file current position to new_posn_B
end

fclose('all');
toc
