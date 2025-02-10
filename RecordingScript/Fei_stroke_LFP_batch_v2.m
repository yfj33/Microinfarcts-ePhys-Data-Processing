clear all
close all
%addpath(genpath('F:\2021-Rice-recording\Fei-UT\2018-stroke\matlab_code'))

% freqrange=[4 30;30 60;60 110];  % this is previous frequency band, no
% indiacation

%Indicate Data recording type
%Recording_Type=128;
Recording_Type=32;

freqrange=[30 60;60 110;300 3000]; % this is new frequency band, indicated by 'new'

window = 1; %unit in seconds.

CMR = 1;
%MinPerFile = 5;  % unit in min
MinPerFile = 1;  % unit in min
HrEachSec = 1; % unit in hr each section
de_Noise = 0;  % denoise on

subFolderList=dir;
subFolderName = {subFolderList.name};
desiredFolder = cellfun(@(y) y(1)=='t' , subFolderName );
subFolderList = subFolderList(desiredFolder);

%%
for n = 1: numel(subFolderList)

    %userpath([pwd '/' subFolderList(n).name]); 
    DIR = dir([pwd '/' subFolderList(n).name '/*.rhd']);
    InitialTime = str2double(DIR(1).name(end-9:end-4));
    str = DIR(1).name(end-16:end-11);
% prompt='InitialTime?: ';
% IntitialTime=input(prompt);
% 
 FilesPerGroup = floor(HrEachSec*60/MinPerFile); % combine each hr's data, HrsEachSec=1,MinPerFile=5;FPG=12
% 
% 
% path1=pwd;
% prompt=([path1([end-4:end-3 end-1:end end-9:end-6]) '?']);
% choice=input(prompt,'s');
% if strcmp(choice,'1')
% str = path1([end-4:end-3 end-1:end end-9:end-6]);
% else
% str = choice  ;  
% end



% %userpath([pwd '/' subFolderList(n).name]); 
cd([pwd '/' subFolderList(n).name]);
for Hr = 1:(numel(DIR)/FilesPerGroup)+1 %DIR gives all the rhd files, numel(DIR) gives number of RHD file,FPG=12
    
for i=(Hr-1)*FilesPerGroup+1:min((Hr)*FilesPerGroup,numel(DIR))
  
  
 read_Intan_RHD2000_file_2021(DIR(i).name);

Fs = frequency_parameters.amplifier_sample_rate;
if i==(Hr-1)*FilesPerGroup+1
data=amplifier_data;
else
    data=[data amplifier_data];
end

% clearvars -except DIR data amplifier_data amplifier_channels frequency_parameters
end


% Do Channel rejection after ch 47, before ch16
y=[amplifier_channels.native_order];imp=[amplifier_channels.electrode_impedance_magnitude];
if Recording_Type==32
selection = (y<48 & y >15)&(imp<2E6); %impedance threshold 2MegaOhms
% if strcmp(str,'10252017')||strcmp(str,'10272017')||strcmp(str,'10292017')
% selection(32:end)=0;
% end
else
    selection=(imp<2E6);
end
amplifier_channels=amplifier_channels(selection);
data=data(selection,:);
%

if de_Noise == 1
    load('denoiseIndex.mat');
    data(:,denoise) = median(median(data));
end

amplifier_data=data; 
window_length = 10000;
dend=floor(size(amplifier_data,2)/window_length)*window_length;
amplifier_data=amplifier_data(:,1:dend);


if CMR ==1
amplifier_data = amplifier_data - repmat(median(amplifier_data),size(amplifier_data,1),1);
elseif CMR==2
    list1=1:15;
    list2=16:31;
amplifier_data(list1,:) = amplifier_data(list1,:) - repmat(median(amplifier_data(list1,:)),size(amplifier_data(list1,:),1),1);
amplifier_data(list2,:) = amplifier_data(list2,:) - repmat(median(amplifier_data(list2,:)),size(amplifier_data(list2,:),1),1);
         
end

%%


 % assume non-overlapping windows. 
 window_length =  floor(window*Fs); %30000 
dend=floor(size(amplifier_data,2)/window_length)*window_length;%dend integer of window number * dp/window
amplifier_data=amplifier_data(:,1:dend); %take part of the data (window number is an integer)
result = cell(size(freqrange,1),1); % result 3*1
% loop through per frequency range per channel
for range=1:numel(result)
   
    perChannelEachWindow = zeros(size(amplifier_data,1),size(amplifier_data,2)/window_length); %Bin dps into window, =channel num* window(1s) num
   for ch=1: size(amplifier_data,1)
   [range ch]
       perChannelEachWindow(ch,:)=bandpower(reshape(amplifier_data(ch,:),window_length,size(amplifier_data,2)/window_length),Fs,freqrange(range,:));%compute avg power in each column
   end%pchanneleachwindow=ch num*binned window#
   result{range}= perChannelEachWindow;
end
x(:,2)=[amplifier_channels.native_order]';
x(:,1)=1:size(x,1);
Dat_V_Map=x;%alive channels 16-47

% movingwin=[5 1];
% params.Fs = Fs;
% params.trialave=0;
% [S,t,f]=mtspecgramc(amplifier_data(8,:)',movingwin,params);
% S=S';
% imagesc(t,f(1:10:3000),S(1:10:3000,:))
% colormap jet 
% caxis([0 10])
% set(gca,'YDir','normal')


% S_normal  = bsxfun(@rdivide,S',median(S(:,1:200)'))';
% imagesc(t,f(1:10:3000),S_normal(1:10:3000,:))
% colormap jet 
% caxis([0 5])
% set(gca,'YDir','normal')

if de_Noise == 1
    save([str '-Hr-' num2str(Hr) 'LFP-new-denoise'],'result','Dat_V_Map','InitialTime');
else
    save([str '-Hr-' num2str(Hr) 'LFP-new'],'result','Dat_V_Map','InitialTime');
end

clear data amplifier_data x
end
cd ..
end
%%
% movingwin=[5 1];
% params.Fs = Fs;
% params.trialave=0;
% [S,t,f]=mtspecgram(amplifier_data(',movingwin,params);
% % save(saveInfo,'S','t','f','params','movingwin')
% cd ..
% % S_low = S(:,1:264,:);
% % S_avg = median(S,3);
% save(saveInfo,'S_low','S_avg', 't','f','params','movingwin','-v7.3')