%% time stamping code
% g3 integration
% Clark Nov to Jan 2020

clear all;
tic

%hist settings
%change these as required- bin size, offset, searching time
t_max=1e6; %time in ps=1ms for searching function
d_bin=100;
bin_size=50*d_bin; %in ps
%t_offset=2.5e6 - bin_size; bunched
t_offset=0;%1.5e6 - bin_size; %antibunched
bins=round(t_max/bin_size);
%full range
histrange=linspace(-t_max,t_max,bins);
N1=zeros(1,length(histrange)-1);
N2=zeros(1,length(histrange)-1);
N3=zeros(length(histrange)-1,length(histrange)-1);
%set up plots
%full range plot
% will have to use surf or contour?
%pick a reference channel - eg C1 - and plot dx=dt(C1-C2), dy=dt(C1-C3)
h=contourf(histrange(1:length(N3)),histrange(1:length(N3)),N3);


%% read in files
%get file parameters such as count rate, number of events
%enter file name here
fileName1='0518_0.9V_100k_10min_C1.bin';
fileName2='0518_0.9V_100k_10min_C2.bin';
fileName3='0518_0.9V_100k_10min_C4.bin';
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
chunk_plot=chunk_factor*1e3;
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
    dt12_array=NaN(ind,200); %preallocate dt for speed
    dt13_array=NaN(ind,200);
    dt_pair=NaN(ind*100,2);
    %dt123_array=NaN(ind,2000);
    i123=0;
    for i2=1:ind
        here=indices(i2); %find index value of this iteration
        now=d_array1(here); %find time value of this index
        hi=now + (2*t_max); %set upper bound
        lo=now - (2*t_max); %set lower bound
        %search in next channel using find
        match2 = find(d_array2<hi & d_array2>lo);
        if isempty(match2) == false
            dt12=d_array2(match2) - now;
            dt12_array(i2,1:length(dt12))=dt12;
        end
        match3 = find(d_array3<hi & d_array3>lo);
        if isempty(match3) == false
            dt13=d_array3(match3) - now;
            dt13_array(i2,1:length(dt13))=dt13;
        end
        
        row12=dt12_array(i2,:);
        row12_t=row12(~isnan(row12));
        i3=length(row12_t);
        row13=dt13_array(i2,:);
        row13_t=row13(~isnan(row13));
        i4=length(row13_t);
        row123=i3*i4;
        
        for i12=1:i3
            for i13=1:i4
                i123=i123+1;
                dt_pair(i123,1)=dt12_array(i2,i12);
                dt_pair(i123,2)=dt13_array(i2,i13);
            %dt123_array(i2,i12*i13)=dt12_array(i2,i12)-dt13_array(i2,i13);
            end
        end
    end
    
    %count data
    [n1, histrange_new]= histcounts(dt12_array,histrange);%create histogram of pairs file- this command is a bottleneck on program speed, non linear with 'edges' size
    N1= N1+n1;
    [n2,histrange_new]=histcounts(dt13_array,histrange);
    N2=N2+n2;
    
    [n3,histrange_new,histrange_new]=histcounts2(dt_pair(:,1),dt_pair(:,2),histrange,histrange);
    for i3=1:length(histrange)-1 %column iterator
        for i4=1:length(histrange)-1 %row iterator
            N3(i4,i3)=N3(i4,i3)+n3(i4,i3);
        end
    end
   
    %% update plot
    i3=mod(i1,round(chunk_plot/chunk_factor));
    if i3==0%only update plot every time chunk_plot_size has been processed- faster
        h.ZDataSource = 'N3';
        refreshdata;
        shg;
    end
    
    T2=isempty(d_array2);
    T3=isempty(d_array3);
    
    if T2 == 1 | T3 == 1
        break
    
    elseif T2 == 0 || T3 == 0
        grad2 = ftell(fid2)/max(d_array2);%gradient of fid2_file_pos at end of this chunk/channel2(time)
        new_posn_2 = 8*round(max(d_array1)*grad2/8);%each timestamp has 8 bytes so find time to nearest 8 units
        fseek(fid2,new_posn_2,'bof');%move file current position to new_posn_B
        grad3=ftell(fid3)/max(d_array3);
        new_posn_3=8*round(max(d_array2)*grad2/8);
        fseek(fid3,new_posn_3,'bof');
    end
    
   
end

% fullfileName=fullfile(fileDir,'N12.dat');
% dlmwrite(fullfileName,N1,'-append');
% fullfileName=fullfile(fileDir,'N13.dat');
% dlmwrite(fullfileName,N2,'-append');
    
fclose('all');
toc


%dtsize(i2)=length(dt1); %size of dt array
    %mini(i2)=min(abs(dt1)); %minimum dt
    %this needs to work such that t12 t13 are same size
    %rather than reinvent hist2, reinvent this calculation step
    
    %remove nans
    %dt_plt12=dt12_array(dt12_array~=0);
    %dt_plt13=dt13_array(dt13_array~=0);
    
    


    %write data to files
%     fileDir= 'C:\Users\rache\Documents\Uni\postgrad\MATLAB';
%     fullfileName=fullfile(fileDir,'dt12.dat');
%     dlmwrite(fullfileName,dt_plt12,'-append');
%     fullfileName=fullfile(fileDir,'dt13.dat');
%     dlmwrite(fullfileName,dt_plt13,'-append');


%     for i3=1:length(histrange)-1 %column iterator
%         for i4=1:length(histrange)-1 %row iterator
%             pos1=N1(i4); %N1 = C1-C2, N2 = C1-C3
%             pos2=N2(i3); %row needs to be C1-C2, column C1-C3
%             if pos1~=0 && pos2~=0
%                 N3(i4,i3)=N1(i4)+N2(i3);
%             end
%         end
%     end
