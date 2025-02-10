%STA Spike Triggered Average 
min_per_file=5; %5min pef file
seconds_a_min=60;
time_window=400; %time window in ms

RHD=dir(pwd);
RHD_name={RHD.name};
rhd_sel=cellfun(@(y) y(1)=='2',RHD_name);
RHD_name_used=RHD_name(rhd_sel);
channel_map=load('pcb_finger_channel_map.mat');
chmap=channel_map.pcb_finger_channel_map;

  n=1;%for n=1:numel(RHD_name_used)
    read_Intan_RHD2000_file_2021(RHD_name_used{n});
    raw_signal=amplifier_data;
    SR=frequency_parameters.amplifier_sample_rate; %sampling rate
    points_num=time_window/1000*SR;
    %CMR
    native_order=cell2mat({amplifier_channels.native_order}); %for channel mapping to location
    Firings=readmda('firings.mda');%Firings event matrix
    
    ClustN=numel(unique(Firings(3,:)));%cluster number
    AllClust=[1:1:ClustN];
    Clust=Firings(3,:);%cluster number
    Event=Firings(2,:);%firing event occurence

    Event_s=Event/SR; % Event in s 
    PriCh=Firings(1,:);%primary channel range from 1-32
    Noise=xlsread('curation.xlsx','A:A');
    Noise=rmmissing(unique(Noise));
    Positive_Spikes=xlsread('curation.xlsx','B:B');
    Positive_Spikes=rmmissing(unique(Positive_Spikes));
    MSUnit=AllClust(~ismember(AllClust,[Noise;Positive_Spikes]));%Multi and Singel units
    
    %Try primary channel spiking events 
    %Fix one reference electrode and compute for all others wave magnitude

    %LFP filter

    %Bandpass 0.3-500 Hz
    signal_band_filtered=bandpass(raw_signal',[0.3 500],SR);
    %Butter worth filtered 3rd high pass 10th low pass 
    %8th order high pass 3Hz
    fhigh=3;
    flow=90;
    [bh,ah]=butter(8,fhigh/(SR/2),"high");
    fvtool(bh,ah);
    signal_hp=filtfilt(bh,ah,signal_band_filtered); %high pass filtered
    %10TH order low pass filter
    [bl,al]=butter(10,flow/(SR/2),"low");
    fvtool(bl,al);
    signal_lp=filtfilt(bl,al,signal_hp);
    %LFP signal z scored
    LFP_preprocessed=zscore(signal_lp); 
    LFP_preprocessed=LFP_preprocessed';
    
    %primary channel spiking extraction
    tstart=(n-1)*min_per_file*seconds_a_min;
    tend=n*min_per_file*seconds_a_min;
    Event_OI=Event_s(tstart<Event_s&Event_s<tend); %Events of interest 
    Event_OI(1)=[]; Event_OI(end)=[];
    PriCh_OI=PriCh(tstart<PriCh&PriCh<tend); %Primary Channel
    PricCh_OI(1)=[];PriCh_OI(end)=[];
    %classify events to their primary channel and compute STA
    Channels=[1:1:unique(PriCh)];
    STA_waveform=zeros(numel(Channels),points_num,numel(Channels));

    for i=1:numel(Channels)
    STA_sum=zeros(1,points_num);
    Event_ChN=Event_OI(PriCh_OI==Channels(i)); %Events at primary channel i in s
    Event_ChN_seq= Event_ChN*SR;

    for s=1:numel(Event_ChN_seq)
       seq= Event_ChN_seq(s);
       sta=LFP_preprocessed(:,seq-(points_num/2-1):seq+point_num/2); %
       STA_sum=sta+STA_sum;
     
    end
   STA_waveform(:,:,i)=STA_sum./numel(Event_ChN_seq);

    end
   
    figure
    plot(STA_waveform(:,:,2));
   
    