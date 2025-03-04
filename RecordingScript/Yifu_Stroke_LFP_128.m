%extract LFPs from target folder and plot them based on recovery week information
%%   Detailed explanation goes here
clc;clear;close all
folder_name='H:\Yifu Jin\LFP Small Scale Stroke\2021-07-20-aged\high fqz Stroke LFPs';
shank_num=[1 2 4]; %compute 5 sets of results with all shanks, shank 1 2 3 4 separately
%save_folder='G:\Yifu\ePhys\2021-09-05-aged\S4';
BL=4;
Chtype=32;

% if ~exist(save_folder,'dir')
%     mkdir(save_folder)
% end
%E:\2021-Rice-recording\2020-21-surgeries\2021-07-28-stroke\Stroke LFPs
%week_info={'-11','-10','-4','Stroke Day','2','3','4','5','6','7','8','10','Week2','Week3','Week4','Week6','Week8'}
%week_info={'-3','-2'}
%week_info={'BL1','BL2','BL3','BL4','Stroke Day','48hrs'};
%week_info={'Baseline 1','Baseline 2','Baseline 3','Baseline 4','Stroke','D2','D3','D4','D5','D6','D7','D10','D14','D21','D28'};
%week_info={'Bl1','Bl2','Bl3','Bl4','S','D2','D3','D4','D5','D6','D7','D10','W2','W3','W4'};
%week_info={'Bl1','Bl2','Bl3','Bl4','S','D2','D3','D4','D5','D6','D7','D10','W2','W3','W4','W5','W6','W7','W8','W9'};
%week_info={'Bl1','Bl2','Bl3','Bl4','Stroke','48hrs','Day5','Day7','Day10','Day14'};

%week_info={'Bl1','Bl2','Bl3','Bl4','S','D2','D3','D4','D5','D6','D7','D8','D10','W2','W3','W4','W6','W8'};
%week_info={'Bl1','Bl2','Bl3','S','D2','D3','D4','D5','D6','D7','D8','W2','W3','W4'};

%freqrange=[0.5 4;4 8;8 12;12 30];  %delta: 0.5-4Hz; theta: 4-8Hz; alpha: 8-12Hz; beta: 12-30Hz
%week_info={'Bl1','Bl2','Bl3','Bl4','SI','D2','D3','D4','D5','D6','D7','D10','W2'};
week_info={'Bl1','Bl2','Bl3','Bl4','SI','D2','D3','D4','D5','D6','D7','D10','W2'};
%week_info={'Bl1','Bl2','Bl3','D1','D2','D3','D4','D5','D6','D7','D10','W2'};
%week_info={'Bl1','Bl2','Bl3','D1','D2','D3','D4','D5','D7','D10','W2'};
if Chtype==32;
    CMN=load('Z:\xl_stroke\yifu_Digigait\SpikeSorting\SpikeSorting_Script\Channel_Map\pcb_finger_channel_map.mat');
    CMN=CMN.pcb_finger_channel_map;
    CMN=CMN';
else
CMN=load('Z:\xl_stroke\yifu_Digigait\SpikeSorting\SpikeSorting_Script\Channel_Map\chan_map_1x32_128ch_rigid.mat');
CMN=CMN.CMN;
CMN=CMN';
end
%CMN=CMN.pcb_finger_channel_map;

%channel_reject=1; %if there is channel to reject put in reject shank/channel number

%if channel_reject==1
reject_ch_num=[CMN(:,[3,1,2])];

%end
    

%% Extratc file name
cd(folder_name)
fileinfo=dir('*.mat');%obtain the file name in current folder
filename={fileinfo.name};%filename list in cells
%% read each LFP file and calculate the average power 
AllData=cell(size(filename,2),1); %prelocate cell to store all the data from all sessions


band_num=3; %3 frequency bands:[30 60;60 110;300 3000]
LFP_avg=nan(size(AllData,1),band_num);
LFP_std=nan(size(AllData,1),band_num); %prelocate: session num*band num for LFP mean and std

LFP_norm_BL=nan(size(AllData,1),band_num);
LFP_std_norm_BL=nan(size(AllData,1),band_num); %prelocate normalized LFP and std
%LFP_norm_L=nan(size(AllData,1),band_num);
%LFP_std_norm_L=nan(size(AllData,1),band_num);%prelocate normalized LFP against last day

%load PCB channel map
%load('F:\2021-Rice-recording\Channel_map\pcb_finger_channel_map.mat');
bin_interval = 30; % unit in second
%co = parula(6); %for plotting

%% loop through each channel and bandnum
%Denoise: get rid of large noise > 75 percentile  +1.5*iqr of
%raw lfp in each channel each band -- Raw lfp


%bin lfp into 30s bin --New LFP
for sk=1:numel(shank_num)+1
for i=1:size(AllData,1)%session num
   for j=1:band_num%
       AllData{i}=load(filename{i});
      
       S_B=AllData{i}.result{j,1}; %go through each band in each session
       %for each session reject shank number 2
       %if channel_reject==1
       %reject_ch_num=pcb_finger_channel_map(3,:);
       if sk==numel(shank_num)+1
        raw_LFP=S_B';
        save_folder=fullfile(folder_name,strcat('All Shanks'));
        sk_info='all shanks';
         else
      % reject_ch_num=[CMN(:,shank_num(s))];
       kept_ch_num=CMN(:,sk); %the channels kept: other 3 shanks
       [C ia ib]=intersect(AllData{i}.Dat_V_Map(:,2),kept_ch_num);
       %noise rejection
       raw_LFP=S_B(ia,:)'; %dps*channel num
       save_folder=fullfile(folder_name,strcat('Shank',num2str(sk)));
        sk_info=num2str(sk);
        end
       if ~exist(save_folder)
        mkdir(save_folder)
      end
       
       
       %noise rejection
       %raw_LFP=S_B(ia,:)'; %dps*channel num
       
       th=prctile(raw_LFP,75)+1.5*iqr(raw_LFP); %threshold set for each channel
       test=bsxfun(@gt,raw_LFP,th);%find the dps greater than threshold,logical 0 and 1
       raw_LFP(test)=median(raw_LFP,'all');%replace noisy point with median of LFP from all channels
       
       %bin noise rejected LFP in 30s and put a in new_LFP matrix
       
      bin_num=floor(size(raw_LFP,1)/bin_interval); %number of bins, integer
      tend=bin_num*bin_interval; %total number of points
      
      new_LFP=squeeze(nanmean(reshape(raw_LFP(1:tend,:),bin_interval,bin_num,size(raw_LFP,2)),1)); %channel num*bin num: average in 30s bin
     
      %LFP average across all channels
      LFP_avg(i,j)=nanmean(new_LFP,'all');%store LFP average across all channels in lfp_avg, one number
      LFP_std(i,j)=std(nanmean(new_LFP,2)); %mean of each channel within 30s and take std across all channels
   end
end

   %normalize against baseline mean value
   for m=1:size(LFP_avg,1)
       for n=1:size(LFP_avg,2)
           LFP_avg_norm_BL(m,n)=LFP_avg(m,n)/mean(LFP_avg(1:BL,n));
           %LFP_avg_norm_BL(m,n)=LFP_avg(m,n)/mean(LFP_avg(1,n));
           LFP_std_norm_BL(m,n)=LFP_std(m,n)/mean(LFP_avg(1:BL,n));
           %LFP_std_norm_BL(m,n)=LFP_std(m,n)/mean(LFP_avg(1,n));
           %LFP_lo_norm_D0(m,n)=LFP_std(m,n)/LFP_avg(1,n);
           %LFP_hi_norm_D0(m,n)=(LFP_avg(m,n)+LFP_std(m,n))/(LFP_avg(1,n)+LFP_std(1,n));
           %normalize against day0
           %LFP_avg_norm_L(m,n)=LFP_avg(m,n)/LFP_avg(end,n);
           %LFP_std_norm_L(m,n)=LFP_std(m,n)/LFP_avg(end,n);
           
           %LFP_lo_norm_L(m,n)=(LFP_avg(m,n)-LFP_std(m,n))/(LFP_avg(end,n)-LFP_std(end,n));
           %LFP_hi_norm_L(m,n)=(LFP_avg(m,n)+LFP_std(m,n))/(LFP_avg(end,n)+LFP_std(end,n));
           
           
       end
       
       
   end
   
  LFP_lo_norm_BL=LFP_avg_norm_BL-LFP_std_norm_BL;
  LFP_hi_norm_BL=LFP_avg_norm_BL+LFP_std_norm_BL;
  %LFP_lo_norm_L= LFP_avg_norm_L-LFP_std_norm_L;
  %LFP_hi_norm_L=LFP_avg_norm_L+LFP_std_norm_L;
  
  
   %get high and low value for D0 and last day normalization
 %LFP_norm_D0_hi=LFP_norm_D0+LFP_std_norm_D0;
 %LFP_norm_D0_lo=LFP_norm_D0-LFP_std_norm_D0;
 
 %LFP_norm_L_hi=LFP_norm_L+LFP_std_norm_L;
 %LFP_norm_L_lo=LFP_norm_L-LFP_std_norm_L;
 
 %save LFP values in the structure
 sp=regexp(folder_name,'\','split');
 mouse_id=sp{end-1};
 SK=struct('LFP_avg_norm_BL',LFP_avg_norm_BL,'LFP_std_norm_BL',LFP_std_norm_BL,'LFP_avg',LFP_avg,'LFP_std',LFP_std,'LFP_lo_norm_BL',LFP_lo_norm_BL,'LFP_hi_norm_BL', LFP_hi_norm_BL);
 %save(strcat(mouse_id,'-LFPs'),'sk');
 
 save(strcat(save_folder,'\',mouse_id,sk_info,'-LFPs'),'SK');%save id in local folder and combined folder
 
 
 color=['r','g','b','y'];
 bandrange={'30 60';'60 110';'300 3000'}%{'0.5-4','4-8','8-12','12-30'};%[30 60;60 110;300 3000] freqrange=[0.5 4;4 8;8 12;12 30;];  %delta: 0.5-4Hz; theta: 4-8Hz; alpha: 8-12Hz; beta: 12-30Hz
%plot of two normalization
%cd(save_folder)
 figure
 wkn=length(week_info);
 x=[1:1:wkn];
 for k=1:band_num
    subplot(band_num,1,k) 
     ss=patch([x fliplr(x)],[ LFP_lo_norm_BL(:,k)' fliplr( LFP_hi_norm_BL(:,k)')],color(k),'EdgeColor',color(k));
 hold on
  plot(1:wkn, LFP_avg_norm_BL(:,k),'color',color(k),'MarkerFaceColor',color(k));
  hold on
  scatter(1:wkn, LFP_avg_norm_BL(:,k),25,color(k),'filled');
  %alpha(s,0.1)
  set(ss,'FaceAlpha',0.25,'edgecolor','none');
  xticklabels(week_info);
  xticks([1:1:numel(week_info)]);
  axis tight
  title(strcat('Baseline Normalized Bandwidth',bandrange{k},'Hz'),'FontSize',28);
  xlabel('Days Relative to Stroke','FontSize',28);
   ylabel('Normalized Value','FontSize',28);
   %ylim([0 max(LFP_avg_norm_BL(:)+0.2)]);
  ax = gca;
ax.XAxis.FontSize = 28;
ax.YAxis.FontSize = 28;
 end
 
%xticklabels(week_info)
%legend('30-60 Hz','60-110 Hz','300-3000Hz')
 %title('LFP normalized against Day0')  
 savefig(fullfile(save_folder,strcat(mouse_id,'Baseline Normalization LFP.fig')));
 
 for k=1:band_num
    figure 
     ss=patch([x fliplr(x)],[ LFP_lo_norm_BL(:,k)' fliplr( LFP_hi_norm_BL(:,k)')],color(k),'EdgeColor',color(k));
 hold on
  plot(1:wkn, LFP_avg_norm_BL(:,k),'color',color(k),'MarkerFaceColor',color(k),'linewidth',1.5);
  hold on
  scatter(1:wkn, LFP_avg_norm_BL(:,k),25,color(k),'filled');
  %alpha(s,0.1)
  set(ss,'FaceAlpha',0.25,'edgecolor','none');
  xticklabels(week_info)
  xticks([1:1:numel(week_info)])
  axis tight
  title(strcat('Baseline Normalized Bandwidth',bandrange{k},'Hz'),'FontSize',28);
  xlabel('Days Relative to Stroke','FontSize',28);
   ylabel('Normalized Value','FontSize',28);
   %ylim([0 max(LFP_avg_norm_BL(:)+0.2)]);
  ax = gca;
ax.XAxis.FontSize = 28;
ax.YAxis.FontSize = 28;
savefig(fullfile(save_folder,strcat(mouse_id,bandrange{k},' Baseline Normalization.fig')));
 end
 close all
end


% figure
 %for k=1:band_num
    % subplot(band_num,1,k)
     
 % v=patch([x fliplr(x)],[ LFP_lo_norm_L(:,k)' fliplr( LFP_hi_norm_L(:,k)')],color(k),'EdgeColor',color(k))
  %hold on
 %  plot(1:wkn,LFP_avg_norm_L(:,k),'color',color(k),'MarkerFaceColor',color(k),'MarkerEdgeColor',color(k))
  % hold on
  % scatter(1:wkn, LFP_avg_norm_L(:,k),25,color(k),'filled')
 %  alpha(v,0.1)
 %  xticklabels(week_info)
 %  title(strcat('Last Day Normalized LFP',bandrange{k}),'FontSize',18)
 % xlabel('Weeks after Surgery','FontSize',18)
 % ylabel('Normalized Value','FontSize',18)
 %  ax = gca;
%ax.XAxis.FontSize = 18;
%ax.YAxis.FontSize = 18;
% end
% savefig(strcat(mouse_id,'Last Day Normalization LFP.fig'))

%xticklabels(week_info)
%legend('30-60 Hz','60-110 Hz','300-3000Hz')
  %title('LFP normalized against Week6')  
    
      
      
      
      




