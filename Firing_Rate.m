%Compute Firing Rate of current animal folder from the metrics information
%Compute cluster number and cluster per channel plot
%same as one section in Manual_Curation
%2022-07-22, compute firing rate on each shank based on primary channel
%2022-08-28 compute temporary single unit firing rate (before merging) over
%final cluster number are subject to change based on merging results
%2022-08-28 add options for animals with severe channel loss: reject channels based on last day remaining channel number
%2022-08-28 it works by counting the channel number on each shank if the
%number drops to more than 50% of its original number then use the channel
%number of last session with (later work)
%First start by manual identify the severe channel loss shank 

%2022-08-28: Output file: FR_by_Shank.mat (FR_Info mean std cluster number),FR NegU per Shank.mat(FR number each each day) 
%FR NegU label per Shank.mat (negative Unit FR labels)  
%Channel Label per Shank.mat (channel number label) 

%2023-01-18 add depth information: FR SU yield amplitude, will need to add
%add multiunits 

clc;clear;close all;
%Options
Channel_Type=32; % for 32 channels type in 32, for 128 type in 128
Ch_Reject_Opt=0; %If 1 compute FR results with channel rejection
Ch_Loss_Shank=[3];% Give the shank to reject manually by looking at the data first

Closest_infarct=[4]; %infarct closest to stroke site less than 600um
%Week information to edit
%week_info={'Bl1','Bl2','Bl3','Stroke','D2','D3','D4','D5','D6','D7','D10','W2','W3','W4','W5','W6'};
BL_num=3;
vioc=2.5; %2%
%week_info={'Bl1','Bl2','Bl3','Bl4','S','D2','D3','D4','D5','D6','D7','D10','W2','W3','W4','W5','W6','W7','W8'};
%week_info={'Bl1','Bl2','Bl3','Bl4','S','D2','D3','D4','D5','D6','D7','D10','W2','W3'};
%week_info={'Bl1','Bl2','Bl3','Bl4','S','D2','D3','D4','D5','D7','D10','W2','W3','W4','W5','W6','W7','W8'};
%week_info={'Bl1','Bl2','Bl3','Bl4','S','D2','D3','D4','D5','D6','D7','D10','W2','W3'}
%week_info={'Bl1','Bl2','Bl3','Bl4','S','D2','D3','D4','D5','D6','D7','D10','W2'};
%x_range=[-4,-2,-1,0,1,2,3,4,5,6,9,13];
x_range=[-4,-3,-2,-1,0,1,2,3,4,5,6,9,13];
% x_range=[-3,-2,-1,0,1,2,3,4,5,6,9,13];
%x_rx_range=[-4,-2,-1,0,1,2,3,4,5,6,9,13];
%x_range=[-4,-3,-2,-1,0,1,2,3,4,5,6,9,13,20,27];

if Channel_Type==32;
    cmap=load('Z:/xl_stroke/yifu_Digigait/SpikeSorting/SpikeSorting_Script/Channel_Map/pcb_finger_channel_map_8x4_32ch.mat');
    CMN=cmap.pcb_finger_channel_map_8x4_32ch;
    
%Load channel map and convert to 1-128
else
    cmap=load('Z:/xl_stroke/yifu_Digigait/SpikeSorting/SpikeSorting_Script/Channel_Map/chan_map_1x32_128ch_rigid.mat');
CMN=cmap.CMN;
end
CMN=CMN+1; %1-128



animal_folder=pwd;
folders=dir;
folders(1:2)=[];
folders_name={folders.name};
desired_folder=cellfun(@(y) y(1)=='2',folders_name);
desired_list=folders_name(desired_folder);
firing_rate_pos=zeros(1,numel(desired_list));
firing_rate_neg=zeros(1,numel(desired_list));
clust=zeros(3,numel(desired_list)); %first row cluster number, second row cluster channel number, third row cluster/channel
%preallocate OUTPUT storage matrix
FR_Neg=cell(numel(desired_list),4); %FR of cluster by 4 shank in a cell
FR_NSU=cell(numel(desired_list),4);%FR of negative single units
FR_Pos=cell(numel(desired_list),4);
ch_num=zeros(numel(desired_list),4);%channel number per shank
FR_Neg_lb=cell(numel(desired_list),4);%to store cluster label
FR_Pos_lb=cell(numel(desired_list),4);
FR_NSU_lb=cell(numel(desired_list),4);%negative single units label by shank
FR_NSU_ch=cell(numel(desired_list),4);%negative single units corresponding primary channel
ch_lb=cell(numel(desired_list),4);%channel label per shank

%depth plot and information for Single Units only
depth_FR_SU=cell(8,4,numel(desired_list));  %FR per channel,or cluster
depth_FR_MS=cell(8,4,numel(desired_list));

depth_amp_p2p=cell(8,4,numel(desired_list));  
depth_amp_absp=cell(8,4,numel(desired_list)); 
depth_SUlb=cell(8,4,numel(desired_list));   
depth_MSlb=cell(8,4,numel(desired_list));
depth_MUlb=cell(8,4,numel(desired_list));

depth_SUlb_num=zeros(8,4,numel(desired_list));
depth_MSlb_num=zeros(8,4,numel(desired_list));
depth_MUlb_num=zeros(8,4,numel(desired_list));
%FR Total of each channel
depth_FR_SU_sum=zeros(8,4,numel(desired_list)); 
depth_FR_MS_sum=zeros(8,4,numel(desired_list)); 
%Unit number Total number per channel
depth_SUnit_sum=zeros(8,4,numel(desired_list)); 
depth_SUFR_weightavg=zeros(8,4,numel(desired_list)); 
depth_MSUnit_sum=zeros(8,4,numel(desired_list)); 
depth_MSUFR_weightavg=zeros(8,4,numel(desired_list)); 


%2ms violation record
vio_rec=cell(1,numel(desired_list));


for i=1:numel(desired_list)
cd(fullfile(pwd,desired_list{i}));
%cluster number
mw=readmda('templates.mda');
unit_num=size(mw,3);
AllClust=[1:1:unit_num];
%curation_result=xlsread('curation.xlsx');
%curation_result=readmatrix('curation.xlsx');
Noise=xlsread('curation.xlsx','A:A');
Noise=rmmissing(unique(Noise));
Positive_Spikes=xlsread('curation.xlsx','B:B');
Positive_Spikes=rmmissing(unique(Positive_Spikes));
% if size(curation_result,2)>1 %if there is positive spikes
% Noise=rmmissing(unique(curation_result(:,1)));
% Positive_Spikes=rmmissing(unique(curation_result(:,2)));
% MSUnit=AllClust(~ismember(AllClust,[Noise;Positive_Spikes]));
% elseif size(curation_result,2)==1 %if there isn't positive spikes
% Noise=rmmissing(unique(curation_result(:,1)));
MSUnit=AllClust(~ismember(AllClust,[Noise;Positive_Spikes]));

clustN=numel(MSUnit); %cluster number: multi and single units
chN=size(mw,1); %channel number 
clust(1,i)=clustN;
clust(2,i)=chN;
clust(3,i)=clustN/chN; %unit/channel

%read sampling rate
%read info JSON for sampling rate
fid_info=fopen('info.json');
sr_info=fread(fid_info,inf);
sr_str=char(sr_info');
fclose(fid_info);
SR=jsondecode(sr_str); %structure with one field
sampling_rate=SR.sample_freq;

%read firing information and sort to shank
Firings=readmda('firings.mda');
Clust=Firings(3,:);%cluster number
Event=Firings(2,:);%firing event occurence
PriCh=Firings(1,:);%primary channel
fid_met=fopen('cluster_metrics.json');
met_info=fread(fid_met,inf);
met_str=char(met_info');
fclose(fid_met);
Metrics=jsondecode(met_str);%read cluster metrics 

%read native channel order and count channel number on each shank
chOrder=readNPY('native_ch_order.npy');
chOrder=chOrder+1; %convert 1-128
%find chOrder within channel map
for k=1:4
    shank_temp=CMN(:,k);
    id=ismember(shank_temp,chOrder);
    ch_lb{i,k}=shank_temp(id);
end

for c=1:numel(chOrder)
[R S]=find(chOrder(c)==CMN);
ch_num(i,S)=ch_num(i,S)+1;
end

%Calculate FR of negative Multi-Single unit clusters and single units alone
%by shank 
%ngn_violation=[];
%SUnit=[];
for m=1:numel(MSUnit)
ngn=MSUnit(m);%negative units label
ngn_pc_og=unique(PriCh(Clust==ngn)); %pc of negative units
ngn_pc=chOrder(ngn_pc_og);
[r,shank]=find(ngn_pc==CMN);
ngn_fr=Metrics.clusters(MSUnit(m)).metrics.firing_rate;
% ngn_amp
ws=mw(ngn_pc_og,:,ngn);
amp_p2p=abs(min(ws)-max(ws));
amp_absp=abs(min(ws));
depth_amp_p2p{r,shank,i}=[depth_amp_p2p{r,shank,i};amp_p2p];
depth_amp_absp{r,shank,i}=[depth_amp_absp{r,shank,i};amp_absp];

ngn_SpikeTrain=Event(Clust==ngn)./sampling_rate*1000; % unit spike train in ms
%ngn_violation=Metrics.clusters(MSUnit(m)).metrics.refractory_violation_2msec;%read violation metrics
ISI_ST=diff(ngn_SpikeTrain); %2ms violation
ngn_violation=numel(ISI_ST(ISI_ST<2))/numel(ISI_ST)*100; %less than 2ms ratio in percentage
vio_rec{1,i}=[vio_rec{1,i} ngn_violation];
depth_FR_MS{r,shank,i}=[depth_FR_MS{r,shank,i};ngn_fr];
depth_MSlb{r,shank,i}=[depth_MSlb{r,shank,i};MSUnit(m)];


if ngn_violation<vioc|ngn_violation==vioc %record single units
    %vio_rec{1,i}=[vio_rec{1,i} ngn_violation];
    FR_NSU_lb{i,shank}=[FR_NSU_lb{i,shank};MSUnit(m)]; %record SU if ISI violation<2%
    FR_NSU{i,shank}=[FR_NSU{i,shank};ngn_fr]; %record SU firing rate 
    FR_NSU_ch{i,shank}=[FR_NSU_ch{i,shank};ngn_pc]; %negative single units primary channel number, multiple units per channel

   depth_FR_SU{r,shank,i}=[depth_FR_SU{r,shank,i};ngn_fr];%firing rate 
   
   depth_SUlb{r, shank,i}=[depth_SUlb{r, shank,i};MSUnit(m)];
    
end   

%FR_Neg=FR_Neg+Metrics.clusters(MSUnit(m)).metrics.firing_rate;
%FR_Neg_1=[FR_Neg_1 Metrics.clusters(MSUnit(m)).metrics.firing_rate];

FR_Neg{i,shank}=[FR_Neg{i,shank};ngn_fr]; %each shank
%FR_Neg_num(i,shank)=FR_Neg_num(i,shank)+1;
FR_Neg_lb{i,shank}=[FR_Neg_lb{i,shank};ngn];
end

cd(animal_folder);
end

%make folder to save information if not exist, folder include structure and
%figures
folder2save=fullfile(animal_folder,'processing results');
if ~exist('folder2save','dir');
    mkdir(folder2save);
end

%calculate and plot depth information
depth_folder=fullfile(folder2save,'depth information');
if ~exist('depth_folder','dir');
    mkdir(depth_folder);
end

for q=1:8
    for e=1:4
        for r=1:numel(desired_list)
           depth_FR_SU_sum(q,e,r)=nansum(depth_FR_SU{q,e,r});
           depth_FR_MS_sum(q,e,r)=nansum(depth_FR_MS{q,e,r});
           %depth_Unit_perCh(q,e,r)=sum(depth_lb{q,:,r}) ;%unit per channel
           depth_SUnit_sum(q,e,r)=numel(depth_SUlb{q,e,r}); %total unit number
           depth_MSUnit_sum(q,e,r)=numel(depth_MSlb{q,e,r});
          % depth_FR_avg(q,e,r)=nanmean(depth_FR{q,e,r});
          weight_matrix_SU=depth_FR_SU{q,e,r}./nansum(depth_FR_SU{q,e,r}) ;
          depth_SUFR_weightavg(q,e,r)= dot(weight_matrix_SU,depth_FR_SU{q,e,r});

          weight_matrix_MS=depth_FR_MS{q,e,r}./nansum(depth_FR_MS{q,e,r}) ;
          depth_MSUFR_weightavg(q,e,r)= dot(weight_matrix_MS,depth_FR_MS{q,e,r});


        end
  
end
end
 % depth_FR_avg(isnan(depth_FR_avg))=0;
%Single Units
SU_depth_info.depth_FR=depth_FR_SU;
SU_depth_info.depth_lb=depth_SUlb;
SU_depth_info.depth_FR_sum=depth_FR_SU_sum;
SU_depth_info.depth_Unit_sum= depth_SUnit_sum;
SU_depth_info.depth_FR_perClu=depth_FR_SU_sum./depth_SUnit_sum; %SU fr PER cluster
SU_depth_info.depth_FR_perCluweight=depth_SUFR_weightavg;
%MultiUnits and Single Units
MSU_depth_info.depth_FR=depth_FR_MS;
MSU_depth_info.depth_lb=depth_MSlb;
MSU_depth_info.depth_FR_sum=depth_FR_MS_sum;
MSU_depth_info.depth_Unit_sum= depth_MSUnit_sum;
MSU_depth_info.depth_FR_perClu=depth_FR_MS_sum./depth_MSUnit_sum;%MS units fr per cluster
MSU_depth_info.depth_FR_perCluweight=depth_MSUFR_weightavg;
MSU_depth_info.depth_amp_p2p=depth_amp_p2p;
MSU_depth_info.depth_amp_absp=depth_amp_absp;
%Multiunits only from previous curated data
depth_FR_MU=cell(8,4,numel(desired_list));
depth_MUlb=cell(8,4,numel(desired_list));
depth_FR_MU_sum=zeros(8,4,numel(desired_list));
depth_MUnit_sum=zeros(8,4,numel(desired_list));

for m=1:8
 for n=1:4
  for i=1:numel(desired_list)
    %multi units fr
   fr_ms=depth_FR_MS{m,n,i};
   fr_su=depth_FR_SU{m,n,i};
   depth_FR_MU{m,n,i}=fr_ms(~ismember(fr_ms,fr_su));
   %multi units label
   lb_ms=depth_MSlb{m,n,i};
   lb_su=depth_SUlb{m,n,i};
   depth_MUlb{m,n,i}=lb_ms(~ismember(lb_ms,lb_su));
   depth_FR_MU_sum(m,n,i)=nansum(depth_FR_MU{m,n,i});
   depth_MUnit_sum(m,n,i)=numel(depth_FR_MU{m,n,i});

  end
 end
end
depth_MUFR_perClu=depth_FR_MU_sum./depth_MUnit_sum;
%MultiUnits only
MU_depth_info.depth_FR=depth_FR_MU;
MU_depth_info.depth_lb=depth_MUlb;
MU_depth_info.depth_Unit_sum=depth_MUnit_sum;
MU_depth_info.depth_FR_sum=depth_FR_MU_sum;
MU_depth_info.depth_FR_perClu=depth_MUFR_perClu;

%save depth information
save(fullfile(depth_folder,'Single Units depth_info.mat'),'SU_depth_info','-mat');
save(fullfile(depth_folder,'MS Units depth_info.mat'),'MSU_depth_info','-mat');
save(fullfile(depth_folder,'MU Units depth_info.mat'),'MU_depth_info','-mat');
%plot 3d heat map FOR single units only
colorlim="auto";
f81=figure %FR per cluster
FR_Xdata=x_range;
FR_Ydata=[0 700];
FR_Cdata= squeeze(nansum(depth_FR_SU_sum,2))./squeeze(nansum(depth_SUnit_sum,2));
FR_Cdata(isnan(FR_Cdata))=0;
%imagesc(FR_Xdata,FR_Ydata,FR_Cdata,'Interpolation','Bilinear');
imagesc(FR_Xdata,FR_Ydata,FR_Cdata)
colormap ("jet");
colorbar;
caxis(colorlim);
title('Firing Rate(Hz)');
ylim([0 800]);
xlabel('Days Relative to Stroke Induction');
ylabel('Depth(µm)');
saveas(f81,fullfile(depth_folder,'Firing Rate at Depth.png'));

f811=figure %FR per cluster for closest infarct 
FR_Xdata=x_range;
%FR_Ydata=[0:100:700];
FRavg_Cdata= squeeze(nansum(depth_FR_SU_sum(:,Closest_infarct,:),2))./squeeze(nansum(depth_SUnit_sum(:,Closest_infarct,:),2));
FRavg_Cdata(isnan(FRavg_Cdata))=0;
imagesc(FR_Xdata,FR_Ydata,FRavg_Cdata);
colormap("jet");
colorbar;
caxis(colorlim);
title('Closest Shank Firing Rate(Hz)');
%ylim([0 800]);
xlabel('Days Relative to Stroke Induction');
ylabel('Depth(µm)');
saveas(f811,fullfile(depth_folder,'Closest Shank Firing Rate at Depth.png'));


f82=figure %unit per channel
U_Xdata=x_range;
U_Ydata=[0:700];
U_Cdata= squeeze(nansum(depth_SUnit_sum,2));
imagesc(U_Xdata,U_Ydata,U_Cdata);
colormap("jet");
colorbar;
caxis(colorlim);
title('Unit/Ch');
%ylim([0 800]);
xlabel('Days Relative to Stroke Induction');
ylabel('Depth(µm)');
saveas(f82,fullfile(depth_folder,'Unit per Ch at Depth.png'));

f821=figure %FR per cluster for closest infarct 
FR_Xdata=x_range;
%FR_Ydata=[0:100:700];
Uavg_Cdata= squeeze(nansum(depth_SUnit_sum(:,Closest_infarct,:),2));
imagesc(FR_Xdata,FR_Ydata,Uavg_Cdata);
colormap("jet");
colorbar;
caxis(colorlim);
title('Closest Shank Unit per Channel');
%ylim([0 800]);
xlabel('Days Relative to Stroke Induction');
ylabel('Depth(µm)');
saveas(f821,fullfile(depth_folder,'Closest Shank Unit per Channel at Depth.png'));

% f81=figure %FR per cluster
% FR_Xdata=x_range;
% FR_Ydata=[0 700];
% FR_Cdata= squeeze(nansum(depth_FR_SU_sum,2))./squeeze(nansum(depth_SUnit_sum,2));
% FR_Cdata(isnan(FR_Cdata))=0;
% imagesc(FR_Xdata,FR_Ydata,FR_Cdata);
% colormap("jet");
% colorbar;
% caxis(colorlim);
% title('Firing Rate(Hz)');
% ylim([0 800]);
% xlabel('Days Relative to Stroke Induction');
% ylabel('Depth(µm)');
% saveas(f81,fullfile(depth_folder,'Firing Rate at Depth.png'));
% 
% f811=figure %FR per cluster for closest infarct 
% FR_Xdata=x_range;
% %FR_Ydata=[0:100:700];
% FRavg_Cdata= squeeze(nansum(depth_FR_SU_sum(:,Closest_infarct,:),2))./squeeze(nansum(depth_SUnit_sum(:,Closest_infarct,:),2));
% FRavg_Cdata(isnan(FRavg_Cdata))=0;
% imagesc(FR_Xdata,FR_Ydata,FRavg_Cdata,'Interpolation','Bilinear');
% colormap("jet");
% colorbar;
% caxis(colorlim);
% title('Closest Shank Firing Rate(Hz)');
% %ylim([0 800]);
% xlabel('Days Relative to Stroke Induction');
% ylabel('Depth(µm)');
% saveas(f811,fullfile(depth_folder,'Closest Shank Firing Rate at Depth.png'));

%
f91=figure %FR per cluster
FR_Xdata=x_range;
FR_Ydata=[0 700];
FR_Cdata= squeeze(nansum(depth_FR_MS_sum,2))./squeeze(nansum(depth_MSUnit_sum,2));
FR_Cdata(isnan(FR_Cdata))=0;
imagesc(FR_Xdata,FR_Ydata,FR_Cdata);
colormap("jet");
colorbar;
caxis(colorlim);
title('Firing Rate(Hz)');
ylim([0 800]);
xlabel('Days Relative to Stroke Induction');
ylabel('Depth(µm)');
saveas(f91,fullfile(depth_folder,'Firing Rate at Depth.png'));

f911=figure %FR per cluster for closest infarct 
FR_Xdata=x_range;
%FR_Ydata=[0:100:700];
FRavg_Cdata= squeeze(nansum(depth_FR_MS_sum(:,Closest_infarct,:),2))./squeeze(nansum(depth_MSUnit_sum(:,Closest_infarct,:),2));
FRavg_Cdata(isnan(FRavg_Cdata))=0;
imagesc(FR_Xdata,FR_Ydata,FRavg_Cdata);
colormap("jet");
colorbar;
caxis(colorlim);
title('Closest Shank Firing Rate(Hz)');
%ylim([0 800]);
xlabel('Days Relative to Stroke Induction');
ylabel('Depth(µm)');
saveas(f911,fullfile(depth_folder,'Closest Shank Firing Rate at Depth.png'));

f92=figure %unit per channel
U_Xdata=x_range;
U_Ydata=[0:700];
U_Cdata= squeeze(nansum(depth_MSUnit_sum,2));
imagesc(U_Xdata,U_Ydata,U_Cdata);
colormap("jet");
colorbar;
caxis(colorlim);
title('Unit/Ch');
%ylim([0 800]);
xlabel('Days Relative to Stroke Induction');
ylabel('Depth(µm)');
saveas(f92,fullfile(depth_folder,'Unit per Ch at Depth.png'));

f921=figure %FR per cluster for closest infarct 
FR_Xdata=x_range;
%FR_Ydata=[0:100:700];
Uavg_Cdata= squeeze(nansum(depth_MSUnit_sum(:,Closest_infarct,:),2));
imagesc(FR_Xdata,FR_Ydata,Uavg_Cdata);
colormap("jet");
colorbar;
caxis(colorlim);
title('Closest Shank Unit per Channel');
%ylim([0 800]);
xlabel('Days Relative to Stroke Induction');
ylabel('Depth(µm)');
saveas(f921,fullfile(depth_folder,'Closest Shank Unit per Channel at Depth.png'));



%save('Firing Rate Positive Units by Shanks','FR_Pos','-mat');
%save('Firing Rate Negative Units by Shanks','FR_Neg','-mat');
%save('Cluster information','clust','-mat');

%plot FR over sessions by shank, preallocate matrix to calculate mean,
%std,cluster count
%FR_Pos_Info=zeros(numel(desired_list),4,3);
FR_Neg_Info=zeros(numel(desired_list),4,3);
FR_NSU_Info=zeros(numel(desired_list),4,3);
for i=1:numel(desired_list)
    for j=1:4
%FR_Pos_Info(i,j,1)=sum(FR_Pos{i,j}); %first z layer for mean value
FR_Neg_Info(i,j,1)=sum(FR_Neg{i,j}); 
FR_NSU_Info(i,j,1)=sum(FR_NSU{i,j});

%FR_Pos_Info(i,j,2)=std(FR_Pos{i,j}); %second z layer for std value
FR_Neg_Info(i,j,2)=std(FR_Neg{i,j});
FR_NSU_Info(i,j,2)=std(FR_NSU{i,j});

%FR_Pos_Info(i,j,3)=numel(FR_Pos{i,j}); %third z layer for cluster count
FR_Neg_Info(i,j,3)=numel(FR_Neg{i,j}); 
FR_NSU_Info(i,j,3)=numel(FR_NSU{i,j});

    end
end
%calculate low high dps
%FR_Pos_lh(:,:,1)=FR_Pos_Info(:,:,1)+FR_Pos_Info(:,:,2);%high
%FR_Pos_lh(:,:,2)=FR_Pos_Info(:,:,1)-FR_Pos_Info(:,:,2);%low
FR_Neg_lh(:,:,1)=FR_Neg_Info(:,:,1)+FR_Neg_Info(:,:,2);%high
FR_Neg_lh(:,:,2)=FR_Neg_Info(:,:,1)-FR_Neg_Info(:,:,2);%low
FR_NSU_lh(:,:,1)=FR_NSU_Info(:,:,1)+FR_NSU_Info(:,:,2);%high
FR_NSU_lh(:,:,2)=FR_NSU_Info(:,:,1)-FR_NSU_Info(:,:,2);%low

%Pos_Folder=fullfile(folder2save,'Positive Spikes by Shank');
Neg_Folder=fullfile(folder2save,'Negative Spikes by Shank');
NSU_Folder=fullfile(folder2save,'Negative Single Units by Shank');

%if ~exist(Pos_Folder,'dir');
%    mkdir(Pos_Folder);
%end

if ~exist(Neg_Folder,'dir');
    mkdir(Neg_Folder);
end

if ~exist(NSU_Folder,'dir');
    mkdir(NSU_Folder);
end



%1st figure: FR per shank for positive units
for s=1:4 %4shank by column
%f1=figure;
%S=patch([x_range fliplr(x_range)],[FR_Pos_lh(:,s,2)' fliplr(FR_Pos_lh(:,s,1)')],'c','EdgeColor','c');
%hold on
%plot(x_range,FR_Pos_Info(:,s,1),'LineWidth',3,'Color','c');
%scatter(x_range,FR_Pos_Info(:,s,1),25,'c','filled');
%xticks(x_range);
%xticklabels(week_info);
%ax=gca;
%ax.XAxis.FontSize = 18;
%ax.YAxis.FontSize = 18;
%set(S,'FaceAlpha',0.25,'edgecolor','none');
%xlabel('Days','FontSize',20);
%ylabel('Firing Rate (Hz)','FontSize',20);
%title(strcat('Positive Units Firing Rate of Shank ',num2str(s)));
%saveas(f1,fullfile(Pos_Folder,strcat('Positive Firing Rate of Shank ',num2str(s),'.fig')));
%close all;

% %f11=figure;%FR per cluster
% S=patch([x_range fliplr(x_range)],[[FR_Pos_lh(:,s,1)./FR_Pos_Info(:,s,3)]' fliplr([FR_Pos_lh(:,s,2)./FR_Pos_Info(:,s,3)]')],'c','EdgeColor','c');
% hold on
% plot(x_range,FR_Pos_Info(:,s,1)./FR_Pos_Info(:,s,3),'LineWidth',3,'Color','c');
% scatter(x_range,FR_Pos_Info(:,s,1)/FR_Pos_Info(:,s,3),25,'c','filled');
% xticks(x_range);
% %xticklabels(week_info);
% ax=gca;
% ax.XAxis.FontSize = 18;
% ax.YAxis.FontSize = 18;
% set(S,'FaceAlpha',0.25,'edgecolor','none');
% xlabel('Days','FontSize',20);
% ylabel('Firing Rate (Hz)','FontSize',20);
% title(strcat('Positive Units Firing Rate per Cluster of Shank ',num2str(s)));
% saveas(f11,fullfile(Pos_Folder,strcat('Positive Firing Rate per Cluster of Shank ',num2str(s),'.fig')));
% close all;
% 
% f12=figure;%FR per channel
% S=patch([x_range fliplr(x_range)],[[FR_Pos_lh(:,s,2)./ch_num(:,s)]' fliplr([FR_Pos_lh(:,s,1)./ch_num(:,s)]')],'c','EdgeColor','c');
% hold on
% plot(x_range,FR_Pos_Info(:,s,1)./ch_num(:,s),'LineWidth',3,'Color','c');
% scatter(x_range,FR_Pos_Info(:,s,1)./ch_num(:,s),25,'c','filled');
% xticks(x_range);
% %xticklabels(week_info);
% ax=gca;
% ax.XAxis.FontSize = 18;
% ax.YAxis.FontSize = 18;
% set(S,'FaceAlpha',0.25,'edgecolor','none');
% xlabel('Days','FontSize',20);
% ylabel('Firing Rate (Hz)','FontSize',20);
% title(strcat('Positive Units Firing Rate per Channel of Shank ',num2str(s)));
% saveas(f12,fullfile(Pos_Folder,strcat('Positive Firing Rate per Channel of Shank ',num2str(s),'.fig')));
% close all;


f2=figure; %negative units
n2=mean(FR_Neg_Info(1:BL_num,s,1));
S=patch([x_range fliplr(x_range)],[FR_Neg_lh(:,s,2)'/n2 fliplr(FR_Neg_lh(:,s,1)'/n2)],'b','EdgeColor','b');
hold on
plot(x_range,FR_Neg_Info(:,s,1)/n2,'LineWidth',3,'Color','b');
scatter(x_range,FR_Neg_Info(:,s,1)/n2,25,'blue','filled');
xticks(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
set(S,'FaceAlpha',0.25,'edgecolor','none');
xlabel('Days','FontSize',20);
ylabel('Firing Rate (Hz)','FontSize',24);
title(strcat('Negative Units Firing Rate of Shank ',num2str(s)));
saveas(f2,fullfile(Neg_Folder,strcat('Negative Firing Rate by Shank ',num2str(s),'.fig')));
close all;

f21=figure;%FR per cluster
n21=mean(FR_Neg_lh(1:BL_num,s,2)./FR_Neg_Info(1:BL_num,s,3));
S=patch([x_range fliplr(x_range)],[[FR_Neg_lh(:,s,2)./FR_Neg_Info(:,s,3)/n21]' fliplr([FR_Neg_lh(:,s,1)./FR_Neg_Info(:,s,3)/n21]')],'b','EdgeColor','b');
hold on
plot(x_range,FR_Neg_Info(:,s,1)./FR_Neg_Info(:,s,3)/n21,'LineWidth',3,'Color','b');
scatter(x_range,FR_Neg_Info(:,s,1)./FR_Neg_Info(:,s,3)/n21,25,'blue','filled');
xticks(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
set(S,'FaceAlpha',0.25,'edgecolor','none');
xlabel('Days','FontSize',24);
ylabel('Firing Rate (Hz)','FontSize',24);
title(strcat('Negative Units Firing Rate per Cluster of Shank ',num2str(s)));
saveas(f21,fullfile(Neg_Folder,strcat('Negative Firing Rate per Cluster of Shank ',num2str(s),'.fig')));
close all;

f22=figure;%FR per channel
n22=mean(FR_Neg_lh(1:BL_num,s,2)./ch_num(1:BL_num,s));
S=patch([x_range fliplr(x_range)],[[FR_Neg_lh(:,s,2)./ch_num(:,s)/n22]' fliplr([FR_Neg_lh(:,s,1)./ch_num(:,s)/n22]')],'b','EdgeColor','b');
hold on
plot(x_range,FR_Neg_Info(:,s,1)./ch_num(:,s)/n22,'LineWidth',3,'Color','b');
scatter(x_range,FR_Neg_Info(:,s,1)./ch_num(:,s)/n22,25,'b','filled');
xticks(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
set(S,'FaceAlpha',0.25,'edgecolor','none');
xlabel('Days','FontSize',24);
ylabel('Firing Rate (Hz)','FontSize',24);
title(strcat('Negative Units Firing Rate per Channel of Shank ',num2str(s)));
saveas(f22,fullfile(Neg_Folder,strcat('Negative Firing Rate per Channel of Shank ',num2str(s),'.fig')));

f6=figure; %negative single units
n6=mean(FR_NSU_Info(1:BL_num,s,1));
S=patch([x_range fliplr(x_range)],[FR_NSU_lh(:,s,2)'/n6 fliplr(FR_NSU_lh(:,s,1)'/n6)],'b','EdgeColor','b');
hold on
plot(x_range,FR_NSU_Info(:,s,1)/n6,'LineWidth',3,'Color','b');
scatter(x_range,FR_NSU_Info(:,s,1)/n6,25,'blue','filled');
xticks(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
set(S,'FaceAlpha',0.25,'edgecolor','none');
xlabel('Days','FontSize',24);
ylabel('Firing Rate (Hz)','FontSize',24);
title(strcat('Negative Single Units Firing Rate of Shank ',num2str(s)));
saveas(f6,fullfile(NSU_Folder,strcat('Negative Single Firing Rate by Shank ',num2str(s),'.fig')));
close all;

f61=figure;%FR per cluster
n61=mean(FR_NSU_Info(1:BL_num,s,1)./FR_NSU_Info(1:BL_num,s,3));
S=patch([x_range fliplr(x_range)],[[FR_NSU_lh(:,s,2)./FR_NSU_Info(:,s,3)/n61]' fliplr([FR_NSU_lh(:,s,1)./FR_NSU_Info(:,s,3)/n61]')],'b','EdgeColor','b');
hold on
plot(x_range,FR_NSU_Info(:,s,1)./FR_NSU_Info(:,s,3)/n61,'LineWidth',3,'Color','b');
scatter(x_range,FR_NSU_Info(:,s,1)./FR_NSU_Info(:,s,3)/n61,25,'blue','filled');
xticks(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
set(S,'FaceAlpha',0.25,'edgecolor','none');
xlabel('Days','FontSize',24);
ylabel('Firing Rate (Hz)','FontSize',24);
title(strcat('Negative Single Units Firing Rate per Cluster of Shank ',num2str(s)));
saveas(f61,fullfile(NSU_Folder,strcat('Negative Single Units Firing Rate per Cluster of Shank ',num2str(s),'.fig')));
close all;

f62=figure;%FR per channel
n62=mean(FR_NSU_Info(1:BL_num,s,1)./ch_num(1:BL_num,s));
S=patch([x_range fliplr(x_range)],[[FR_NSU_lh(:,s,2)./ch_num(:,s)/n62]' fliplr([FR_NSU_lh(:,s,1)./ch_num(:,s)/n62]')],'b','EdgeColor','b');
hold on
plot(x_range,FR_NSU_Info(:,s,1)./ch_num(:,s)/n62,'LineWidth',3,'Color','b');
scatter(x_range,FR_NSU_Info(:,s,1)./ch_num(:,s)/n62,25,'b','filled');
xticks(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
set(S,'FaceAlpha',0.25,'edgecolor','none');
xlabel('Days','FontSize',24);
ylabel('Firing Rate (Hz)','FontSize',24);
title(strcat('Negative Single Units Firing Rate per Channel of Shank ',num2str(s)));
saveas(f62,fullfile(NSU_Folder,strcat('Negative Single Units Firing Rate per Channel of Shank ',num2str(s),'.fig')));

close all;


end

%f2=figure;
%plot(x_range,firing_rate_neg,'LineWidth',8,'Color','b');
%xticks(x_range);
%xticklabels(week_info);
%ax=gca;
%ax.XAxis.FontSize = 16;
%ax.YAxis.FontSize = 16;
%xlabel('Days','FontSize',24);
%ylabel('Firing Rate (Hz)','FontSize',24);
%saveas(f2,'Firing Rate Plot of Negative Units');



%plot number of units and channel per day 
f3=figure;
plot(x_range,clust(1,:),'LineWidth',3,'Color','b');
hold on
plot(x_range,clust(2,:),'LineWidth',3,'Color','r');
%plot(x_range,clust(3,:),'LineWidth',10,'Color','g');
legend('Cluster Number','Channel Number');
xticks(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
xlabel('Days','FontSize',24);
ylabel('Count','FontSize',24);
saveas(f3,fullfile(folder2save,'Cluster Number.fig'));

%plot cluster per channel 
f4=figure;
plot(x_range,clust(3,:),'LineWidth',8,'Color','g');
%plot(x_range,clust(3,:),'LineWidth',10,'Color','g');
%legend('Cluster Number','Channel Number','Clusters per Channel');
xticks(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
xlabel('Days','FontSize',24);
ylabel('Count','FontSize',24);
saveas(f4,fullfile(folder2save,'Cluster Number per channel.fig'));

%f5  firing rate per units and channels for firing rate of all 4 shanks

FR_Neg_sum=zeros(numel(desired_list),1); %sum of all 4 shanks
%FR_Pos_sum=zeros(numel(desired_list),1);

for i=1:numel(desired_list)
FR_Neg_sum(i,1)=sum(cell2mat(FR_Neg(i,:)'));
%FR_Pos_sum(i,i)=sum(cell2mat(FR_Pos(i,:)'));
end

% f51=figure;
% plot(x_range,FR_Pos_sum);
% xticks(x_range);
% ax=gca;
% ax.XAxis.FontSize = 18;
% ax.YAxis.FontSize = 18;
% xlabel('Days','FontSize',20);
% ylabel('FR(Hz)','FontSize',20);
% title('Positive Units Firing Rate');
% saveas(f51,fullfile(folder2save,'Positive Units Firing Rate.fig'));

f52=figure;
plot(x_range,FR_Neg_sum,'Linewidth',8);
xticks(x_range);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
xlabel('Days','FontSize',24);
ylabel('FR(Hz)','FontSize',24);
title('Negative Units Firing Rate');
saveas(f52,fullfile(folder2save,'Negative Units Firing Rate.fig'));

f53=figure;
plot(x_range,FR_Neg_sum./clust(1,:)','Linewidth',8);
xticks(x_range);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize =24;
xlabel('Days','FontSize',24);
ylabel('FR(Hz)','FontSize',24);
title('Negative Units Firing Rate per Cluster');
saveas(f53,fullfile(folder2save,'Negative Units Firing Rate per Cluster.fig'));

%Single Units Plot by shank

%save FR of each shank:mean std cluster count and channel number of each
%day
%FR_by_Shank=struct('FR_Pos',FR_Pos,'FR_Pos_Info',FR_Pos_Info,'FR_Neg',FR_Neg,'FR_Neg_Info',FR_Neg_Info);
FR_by_Shank_info=struct('FR_Neg_Info',FR_Neg_Info,'FR_NSU_Info',FR_NSU_Info,'ch_num',ch_num,'ch_lb',ch_lb);
save(fullfile(folder2save,'FR_by_Shank.mat'),'FR_by_Shank_info');
save(fullfile(folder2save,"Cluster and Channels"),'clust');
%save cluster firing, cluster label, repective FR ,channel label of each shank in separate
%files
save(fullfile(Neg_Folder,'FR NegU per Shank.mat'),'FR_Neg');
save(fullfile(Neg_Folder,'FR NegU label per Shank.mat'),'FR_Neg_lb');
%save(fullfile(folder2save,'Channel Label per Shank.mat'),'ch_lb');
save(fullfile(NSU_Folder,'FR SingleUnit per Shank.mat'),'FR_NSU');
save(fullfile(NSU_Folder,'FR SingleUnit label per Shank.mat'),'FR_NSU_lb');
save(fullfile(NSU_Folder,'FR Singleunit channel label per Shank.mat'),'FR_NSU_ch');


%find common channels for heavy channel loss shanks and reject some
%clusters or if all shanks have heavy channel loss at someday discard days
%after that or enter all 4 shanks
if Ch_Reject_Opt==1
NSU_Loss_Folder=fullfile(folder2save,'Negative Single Units by Shank with Channel Loss');  

if ~exist(NSU_Loss_Folder)
    mkdir(NSU_Loss_Folder);
end
FR_NSU_Loss=cell(numel(desired_list),4);%FR of negative single units for channel loss
FR_NSU_Loss_clu=cell(numel(desired_list),4);%Cluster label of negative single units for channel loss
FR_NSU_Loss_ch=cell(numel(desired_list),4);%Channel label of negative single units
FR_NSU_Loss_Info=zeros(numel(desired_list),4,3); %info of mean std and cluster number


All_Shank=[1 2 3 4];
Good_Shank=All_Shank(~ismember(All_Shank,Ch_Loss_Shank));
ch_num_loss=ch_num(:,Ch_Loss_Shank);

FR_NSU_Loss(:,Good_Shank)=FR_NSU(:,Good_Shank);
FR_NSU_Loss_clu(:,Good_Shank)=FR_NSU_lb(:,Good_Shank);
FR_NSU_Loss_ch(:,Good_Shank)=FR_NSU_Loss_ch(:,Good_Shank);


for S=1:numel(Ch_Loss_Shank)
S_loss=Ch_Loss_Shank(S);
ch_temp=ch_num(:,S_loss);
min_nonzero=min(ch_temp(ch_temp>0));
min_id=find(ch_temp==min_nonzero);%which day has the most loss
SU_chlb=FR_NSU_ch{min_id,S_loss};   %find the channel lbs from SU at that day and use it as common channels across days
for K=1:numel(desired_list)
    allSU_ch=FR_NSU_ch{K,S_loss};
    [cmmn,allsu_id,part_id]=intersect(allSU_ch,SU_chlb);
    FR_temp=FR_NSU{K,S_loss};
    FR_NSU_Loss{K,S_loss}=FR_temp(allsu_id);

    Clu_temp=FR_NSU_lb{K,S_loss};
    FR_NSU_Loss_clu{K,S_loss}=Clu_temp(allsu_id);

    Ch_temp=FR_NSU_ch{K,S_loss};
    FR_NSU_Loss_ch{K,S_loss}= Ch_temp(allsu_id);
end



end
for i=1:numel(desired_list)
    for j=1:4

FR_NSU_Loss_Info(i,j,1)=sum(FR_NSU_Loss{i,j});%mean



FR_NSU_Loss_Info(i,j,2)=std(FR_NSU_Loss{i,j});%std


FR_NSU_Loss_Info(i,j,3)=numel(FR_NSU_Loss{i,j});%units number

    end
end
%calculate low high dps
FR_NSU_Loss_lh(:,:,1)=FR_NSU_Loss_Info(:,:,1)+FR_NSU_Loss_Info(:,:,2);%high
FR_NSU_Loss_lh(:,:,2)=FR_NSU_Loss_Info(:,:,1)-FR_NSU_Loss_Info(:,:,2);%low

%Plot 
for s=1:4
f7=figure; %negative single units
n7=mean(FR_NSU_Loss_Info(1:BL_num,s,1));
S=patch([x_range fliplr(x_range)],[FR_NSU_Loss_lh(:,s,2)'/n7 fliplr(FR_NSU_Loss_lh(:,s,1)')/n7],'b','EdgeColor','b');
hold on
plot(x_range,FR_NSU_Loss_Info(:,s,1)/n7,'LineWidth',3,'Color','b');
scatter(x_range,FR_NSU_Loss_Info(:,s,1)/n7,25,'blue','filled');
xticks(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
set(S,'FaceAlpha',0.25,'edgecolor','none');
xlabel('Days','FontSize',24);
ylabel('Firing Rate (Hz)','FontSize',24);
title(strcat('Negative Single Units with Channel Loss Firing Rate of Shank ',num2str(s)));
saveas(f7,fullfile(NSU_Loss_Folder,strcat('Negative Single Firing Rate by Shank ',num2str(s),'.fig')));
close all;

f71=figure;%FR per cluster
n71=mean(FR_NSU_Loss_Info(1:BL_num,s,1)./FR_NSU_Loss_Info(1:BL_num,s,3));
S=patch([x_range fliplr(x_range)],[[FR_NSU_Loss_lh(:,s,2)./FR_NSU_Loss_Info(:,s,3)/n71]' fliplr([FR_NSU_Loss_lh(:,s,1)./FR_NSU_Loss_Info(:,s,3)/n71]')],'b','EdgeColor','b');
hold on
plot(x_range,FR_NSU_Loss_Info(:,s,1)./FR_NSU_Loss_Info(:,s,3)/n71,'LineWidth',3,'Color','b');
scatter(x_range,FR_NSU_Loss_Info(:,s,1)./FR_NSU_Loss_Info(:,s,3)/n71,25,'blue','filled');
xticks(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
set(S,'FaceAlpha',0.25,'edgecolor','none');
xlabel('Days','FontSize',24);
ylabel('Firing Rate (Hz)','FontSize',24);
title(strcat('Negative Single Units with Channel Loss Firing Rate per Cluster of Shank ',num2str(s)));
saveas(f71,fullfile(NSU_Loss_Folder,strcat('Negative Single Units with Channel Loss Firing Rate per Cluster of Shank ',num2str(s),'.fig')));
close all;

f72=figure;%FR per channel
n72=mean(FR_NSU_Loss_Info(1:BL_num,s,1)./ch_num(1:BL_num,s));
S=patch([x_range fliplr(x_range)],[[FR_NSU_Loss_lh(:,s,2)./ch_num(:,s)/n72]' fliplr([FR_NSU_Loss_lh(:,s,1)./ch_num(:,s)/n72]')],'b','EdgeColor','b');
hold on
plot(x_range,FR_NSU_Loss_Info(:,s,1)./ch_num(:,s)/n72,'LineWidth',3,'Color','b');
scatter(x_range,FR_NSU_Loss_Info(:,s,1)./ch_num(:,s)/n72,25,'b','filled');
xticks(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
set(S,'FaceAlpha',0.25,'edgecolor','none');
xlabel('Days','FontSize',24);
ylabel('Firing Rate (Hz)','FontSize',24);
title(strcat('Negative Single Units Firing Rate per Channel of Shank ',num2str(s)));
saveas(f72,fullfile(NSU_Loss_Folder,strcat('Negative Single Units with Channel Loss Firing Rate per Channel of Shank ',num2str(s),'.fig')));
end




save(fullfile(NSU_Loss_Folder,'FR SingleUnit per Shank with Channel Loss.mat'),'FR_NSU_Loss');
save(fullfile(NSU_Loss_Folder,'FR SingleUnit label per Shank with Channel Loss.mat'),'FR_NSU_Loss_clu');
save(fullfile(NSU_Loss_Folder,'FR Singleunit channel label per Shank with Channel Loss.mat'),'FR_NSU_Loss_ch');
save(fullfile(NSU_Loss_Folder,'FR Singleunit Info.mat'),'FR_NSU_Loss_Info')

end
close all;







