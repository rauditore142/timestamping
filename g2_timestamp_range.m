%multiple g2 run
clear all;
tic

%hist settings
t_max=1e8; %time in ps, 1e8ps=1ms for searching function
bin_size=5e2; %in ps
%t_offset=2.5e6 - bin_size; %bunched
t_offset=1.5e6 - bin_size; %antibunched
%t_offset=0; %generated data
bins=round(t_max/bin_size);
%full range
histrange=linspace(-t_max,t_max,bins);
co=zeros(1,length(histrange)-1);

%set up plots
h1=plot(histrange(1:length(co)),co,'.-');
xlabel('Time in ps');
ylabel('Counts');

%% read in files
%get file parameters such as count rate, number of events
%enter file name here
fileName1='g4_850nm_1.2V_0.68_10min_C1.bin';
fileName2='g4_850nm_1.2V_0.68_10min_C2.bin';

% try 68 intensity no LP filter 
%fid1 = fopen(fileName1,'r+');
%fid2=fopen(fileName2,'r+');

t_ap=[1e4,5e4,1e5,5e5,1e6];
trange=length(t_ap);

for j=1:trange
    dt_ap=t_ap(j);
    [fid1,fid2]=remove_after_pulse(fileName1,fileName2,dt_ap,bin_size);
    %channel 1
    fseek(fid1,-8,'eof');%go to end of file
    last1=fread(fid1,8,'uint64');
    events1= ftell(fid1)/8; %total events
    rate1=events1/last1; %channel rate
    t1=(1/rate1); %average spacing of events in channel
    tlast1=last1-t1; %the final event in channel
    %channel 2
    fseek(fid2,-8,'eof');
    last2=fread(fid2,8,'uint64');
    events2=ftell(fid2)/8; %total events
    rate2=events2/last2; %channel rate
    t2=(1/rate2); %average spacing of events in channel
    tlast2=last2-t2; %the final event in channel
    tdif=tlast2-tlast1; %the difference in time between final events in each channel
    
    %initialise arrays
    %specify chunk length to read in to workspace from file
    %5e4 recommended
    chunk_factor=5e4;
    chunk1=(events1/chunk_factor);
    chunk2=(events2/chunk_factor);
    precision = 'uint64';
    chunk_plot=chunk_factor*1e2;
    N_exp=zeros(chunk_factor,1);
    co_exp=zeros(chunk_factor,1);
    t_int=zeros(chunk_factor,1);
    mini=zeros(chunk_factor,1);
    
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
        dt=zeros(ind,100); %preallocate dt for speed
        mini_dt=zeros(ind,1);
        t_int(i1)=(d_array1(end)+d_array2(end))/2;
        %N_exp(i1)=rate1*rate2*binedge*t_int(i1);
        co_exp(i1)=rate1*rate2*bin_size*t_int(i1);
        
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
                mini_dt(i2)=min(abs(difference)); %minimum dt
            end
        end
        mini(i1)=mean(mini_dt);
        %remove nans
        dt_plt=dt(dt~=0);
        
        %count data
        [counts,histrange]=histcounts(dt_plt,histrange);
        co=co+counts;
        
        
        %% update plot
        i3=mod(i1,round(chunk_plot/chunk_factor));
        if i3==0%only update plot every time chunk_plot_size has been processed- faster
            x1=histrange(1:length(co));
            h1.XDataSource = 'x1';
            h1.YDataSource = 'co';
            refreshdata;
            shg;
        end
        
        grad2 = ftell(fid2)/max(d_array2);%gradient of fid2_file_pos at end of this chunk/channel2(time)
        new_posn = 8*round(max(d_array1)*grad2/8);%each timestamp has 8 bytes so find time to nearest 8 units
        fseek(fid2,new_posn,'bof');%move file current position to new_posn_B
    end
    
    t1=-0.5e8;
    x1=find(histrange<t1,1,'last');
    tnorm=mean(co(1:x1));
    co_norm=co/tnorm;
    
    fName=sprintf('0.68_C1C2_%.0g.mat',dt_ap);
    save(fName,'co','x1','co_norm','histrange');
    figName=sprintf('0.68_C1C2_%.0g.fig',dt_ap);
    saveas(gcf,figName);
    
end

fclose('all');
toc