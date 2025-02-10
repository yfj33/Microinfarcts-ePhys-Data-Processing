% the script aims to plot the output of mountainsort after initial Spike
% Sorting
% mountainview GUI function replicate in a more clear presentation 
% average waveform on the all channels, primary channel waveform, ISI histogram %
%Only view channels of adjacent 2 shanks to reduce time 2022-08-13 yj33
clc;clear;close all
cmap=load('Z:\xl_stroke\yifu_Digigait\SpikeSorting\SpikeSorting_Script\Channel_Map\chan_map_1x32_128ch_rigid.mat');
files=dir;
files(1:2)=[];
files_name={files.name};
desired_files=cellfun(@(y) y(1)=='2',files_name);
desired_list=files(desired_files);

for dl=1:numel(desired_list)
cd(fullfile(pwd,desired_list(dl).name));

mkdir(fullfile(pwd,'clusters figures'));
save_path=fullfile(pwd,'c lusters figures');
%%Reorganize the firing event according to cluster number
%files=dir('*.mda');
%fname=files.name;
%readma files for every five minutes and generate average waveform 
%load channel map for electrode designs


%CMN=cmap.Ch_Map_new; %channel map order by shank 
CMN=cmap.CMN;
CMN=CMN+1
CMN_R=reshape(CMN,[],16); %reshape channel to 8*16 so 2 column=1 shank from left to right
%CMN_R=CMN_R+1;
m=8;
n=16;

%read info JSON for sampling rate
fid_info=fopen('info.json');
sr_info=fread(fid_info,inf);
sr_str=char(sr_info');
fclose(fid_info);
SR=jsondecode(sr_str); %structure with one field
sampling_rate=SR.sample_freq;

%read cluster metrics 
fid_met=fopen('combine_metrics_new.json');
met_info=fread(fid_met,inf);
met_str=char(met_info');
fclose(fid_met);
Metrics=jsondecode(met_str);
FN=fieldnames(Metrics.clusters(1).metrics);
numFN=numel(FN);%number of field name


T_Span=30; %unit in min
%load native channel order to fit into channel map
ch_native=readNPY('native_ch_order.npy');
ch_native=ch_native+1; % channel order offset by 1: from 1-128
ch_alive=numel(ch_native);

%load templates as the mean waveform
mw=readmda('templates.mda');




%Data=arrayDatastore(readmda('filt.mda'));
Firings=readmda('firings.mda');%Firings event matrix
ClustN=numel(unique(Firings(3,:)));%cluster number
Clust=Firings(3,:);%cluster number
Event=Firings(2,:);%firing event occurence
PriCh=Firings(1,:);%primary channel

%preallocate cell Fir to store Firing event time stamps at each cluster
Fir=cell(ClustN,1);
ClustList=unique(Firings(3,:));%Cluster List based on sequence small to large
binrange_isi=[0:0.5:50]; %bin size should be 0.5ms to visualize 2ms vilation
binrange_acg=[-50:0.5:50]; 
InT=[1:1:size(mw,2)]./sampling_rate*1000; %in ms
T=T_Span*60;  %min to seconds
binrange_fr=[0:1:T]; %time bin 1s
xlab_isi=[0:10:50];
xlab_acg=[-50:10:50];
parfor i=1:ClustN
    f=figure('Position',get(0,'Screensize'),'Visible','off');
    t=tiledlayout(m+3,n);
    t.TileSpacing='tight';
    t.TileIndexing='rowmajor';
    %t=tiledlayout('flow','TileSpacing','tight');
     Fir{i,1}=Event(Clust==ClustList(i)); %output all the firing events for ith cluster and store in cell array
     PC=unique(PriCh(ClustList(i)==Clust));                            %Find the primary channel number
     PC_map=ch_native(PC) ;          %map relative channel orders in firing to native channel order
    %[h l]=find(PC_map==CMN_R); % column as shank number
    [h l]=find(PC_map==CMN);
    %shank=CMN_R(:,l); %primary channel shank
    shank=CMN(:,l);
    shank_alive=shank(ismember(shank,ch_native)); %pc channel shank which is alive
    a=find(PC_map==ch_native); %find pc and its waveform
    maxV=max(abs(mw(a,:,i)));%maximum value of primary channel
     %plot avg waveform
    %for j=1:ch_alive ;%total available channels
     for j=1:numel(shank_alive)
         
       if shank_alive(j)==PC_map; %plot primary channel number as red 
    [r c]=find(shank_alive(j)==CMN_R);
    seq=(r-1)*n+c; %convert coordinates to index read left to rigtht by rows
    %seq=(r-1)*2*n+2*c;
    %nexttile(seq,[2 2]);
    nexttile(seq);
    %shank_alive(j);
    %subplot(m+3,n,seq);
    p=find(shank_alive(j)==ch_native); %convert back to original order
    plot(InT,mw(p,:,i),'Color','r','LineWidth',1.5);%plot around 3ms interval time average waveform, filled in 32*4 plot
    pbaspect([3 2 1]);
     xlim([0 3.4]);
     
    %ylim([-200 200]);
    ylim([-1.1*maxV 1.1*maxV]);
    ylabel('V (μV)');
    xlabel('T (ms)');
    %title(strcat('Channel ',num2str(ch_native(j)),' waveform'));
    %pbaspect([1 1 1]);
       else
     [r c]=find(shank_alive(j)==CMN_R);
     seq=(r-1)*n+c;
     nexttile(seq);
     %subplot(m+3,n,seq);
      p=find(shank_alive(j)==ch_native); %convert back to original order
     plot(InT,mw(p,:,i),'Color','b','LineWidth',1.5); 
     pbaspect([3 2 1]);
    xlim([0 3.4]);
    %ylim([-200 200]);
     %maxV=max(mw(p,:,i));%maximum value of primary channel
    ylim([-1.1*maxV 1.1*maxV]);
    ylabel('V (μV)');
    xlabel('T (ms)');
       end
    end
    
    %plot ISI/ACG, firing rate plot, 
    ST=Fir{i,1}./sampling_rate.*1000; %convert to ms
    ISI=diff(ST);
    ACG=[ISI,-ISI];%sudo-ACG just for visualization
    
    %ISI Plot
    nexttile(m*n+1,[3 3]);
   % subplot(m+3,n,[m*n+1,3]);
    histogram(ISI,binrange_isi,'FaceColor','b','EdgeColor','none');
    xticks(xlab_isi);
    title('ISI hisograms');
    xlabel('ISI (ms)');
    ylabel('Count');
    axis auto;

    
    %ACG Plot
   nexttile(m*n+5,[3 3]);
    %acg=ACG(ST) %compute acg values
    %subplot(m+3,n,[m*n+5,3]);
    histogram(ACG,binrange_acg,'FaceColor',[0.4940 0.1840 0.5560],'EdgeColor','none'); %plot ACG
    xticks(xlab_acg);
    title('ACG Part');
    xlabel('ACG (ms)');
    ylabel('Count');
    axis auto;
    
    
    %Firing Rate Plot
    nexttile(m*n+9,[3 3]);
    %subplot(m+3,n,[m*n+9,3]);
    st=ST./1000;   %Convert to s
   histogram(st,binrange_fr);
   title('Firing Rate hisograms');
    xlabel('Time (s)');
    ylabel('FR (Hz)');
    axis auto
   
    %Read cluster metrics and display metrics infomation
    nexttile(m*n+13,[3 3]);
    %subplot(m+3,n,[m*n+13,3]);
    CM=Metrics.clusters(i).metrics; %structure with several criteria
    title('Cluster Metrics');
    xlim([0 40]);
    ylim([0 97]);
    for k=1:numFN
       values=getfield(CM,FN{k});
       Words=strcat(FN{k},' : ',num2str(values));
       text(1,2+k*8,Words,'Interpreter','none');
    end
    
  
    saveas(f,fullfile(save_path,strcat('Cluster ',num2str(i))),'png') ;
    %saveas(f,fullfile(save_path,strcat('Cluster ',num2str(i))),'fig') ;
    
close all    
end
clear Firings;
clear mw;
clear PriCh;
clear Metrics;
clear met_info;
clear Metrics;
cd ..
end







% %% read MDA file: 5min section; 
% Min=5;SR=30000;SPM=60; %5min for each processing, Sampling rate=30k, Seconds per min is 60
% NPts=Min*SR*SPM; %Total dps each reading
% F=fopen('filt.mda','rb'); %F is the file ID, 1 std output
% %fread(F,1,dim_type_str); read 4 types of info from headers (int32)
% try
% code=fread(F,1,'int32');
% catch 
%     error('Problem reading file: filt.mda');
% end
% %1st info
% if code>0 
%     num_dims=code;
%     code=-1;
% else
%     fread(F,1,'int32');
%     num_dims=fread(F,1,'int32');
% end
% 
% dim_type_str='int32';
% if num_dims<0
%     num_dims=-num_dims;
%     dim_type_str='int64';
% end
% 
% %preallocate for store dimension array
% S=zeros(1,num_dims);
% for j=1:num_dims
%     S(j)=fread(F,1,dim_type_str); 
% end
% %N=prod(S); %total number of dps
% 
% %preallocate for storage of 5min filtered data
% if num_dims==1 %only one electrode
%     Data=tall(zeros(1,NPts));
% else
%     Data=tall(zeros(S(1),NPts));
% end
% sc=ceil(S(2)/NPts);% 5min sections number
% 
% % Average waveform computation at all available channels using filt.mda,
% % which can be read from templates
% for i=1:sc
% fseek(F,(i-1)*NPts,'bof');%process 5min of data each time, offset at 0, 5min, 10min
% if code==-1 %go through data and read 5min data each time
%     M=zeros(1,NPts*2);
%     M(:)=fread(F,NPts*2,'float');
%     Data(:)=M(1:2:NPts*2)+i*M(2:2:NPts*2); %imaginary part
%  elseif (code==-2)
%     Data(:)=fread(F,NPts,'uchar');
% elseif (code==-3)
%     Data(:)=fread(F,NPts,'float');%code -3 
% elseif (code==-4)
%     Data(:)=fread(F,NPts,'int16');
% elseif (code==-5)
%     Data(:)=fread(F,NPts,'int32');
% elseif (code==-6)
%     Data(:)=fread(F,NPts,'uint16');
% elseif (code==-7)
%     Data(:)=fread(F,NPts,'double');
% elseif (code==-8)
%     Data(:)=fread(F,NPts,'uint32');
% else
%     error('Unsupported data type code: %d',code);
% end;  
% %Average waveform at each firing event, avg 3ms at all available channels
% 
% 
% 
% 
% 
% end












