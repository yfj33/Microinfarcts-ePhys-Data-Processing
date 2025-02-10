clc;clear;close all
%Behavior state: separate moving and static video and compute FR and
%Celltype

%Output:
%For Behavior Video: 1.Moving vs Static firing rate for MSUnits tell whether
%stimulus locked. 2.neuron type during moving and static state based on ACG. 

%make directory to save figures
save_path=fullfile(pwd,'Behavior State Analysis figures');
if ~exist(save_path);
    mkdir(save_path);
end

%read info JSON for sampling rate
fid_info=fopen('info.json');
sr_info=fread(fid_info,inf);
sr_str=char(sr_info');
fclose(fid_info);
SR=jsondecode(sr_str); %structure with one field
sampling_rate=SR.sample_freq;

%Firing rate timestamp
Firings=readmda('firings.mda');%Firings event matrix
ClustN=numel(unique(Firings(3,:)));%cluster number
Clust=Firings(3,:);%cluster number
Event=Firings(2,:);%firing event occurence
PriCh=Firings(1,:);%primary channel
%MS Units Label
%cluster number
mw=readmda('templates.mda');
unit_num=size(mw,3);
AllClust=[1:1:unit_num];
%curation_result=xlsread('curation.xlsx');
curation_result=readmatrix('curation.xlsx');
if size(curation_result,2)>1 %if there is positive spikes
Noise=rmmissing(unique(curation_result(:,1)));
Positive_Spikes=rmmissing(unique(curation_result(:,2)));
MSUnit=AllClust(~ismember(AllClust,[Noise;Positive_Spikes]));
elseif size(curation_result,2)==1 %if there isn't positive spikes
Noise=rmmissing(unique(curation_result(:,1)));
MSUnit=AllClust(~ismember(AllClust,Noise));
end

%Time Stamps
Moving=[0 765];%bound in s
Static=[765 1488];

%Cell type, add path later
cell_type=load('cell_type.mat');
spikes=cell_type.cell_type.spikes;

%%unpack
ts=spikes.ts;
times=spikes.times;
MSUnit=spikes.cluID;
shankID=spikes.shankID;
filtWaveform=spikes.filtWaveform;
spindices=spikes.spindices;
numcells=spikes.numcells;
total=spikes.total;

%%pack into 2 different states
%process 5min data and calculate avg waveform
%[S,A] = readmda_block(10,fname,sr); %10min block

%create matrix to store neuron type  spontaneous excitatory inhibitory as 1
%2 and 3 unders stim
neuron_type_stim=zeros(numel(MSUnit),2);% 1st column unit label, 2,3,4 column neuron type under stim (1 2 3)
neuron_type_stim(:,1)=MSUnit';

for i=1:numel(MSUnit)
    f=figure('visible','off');
    box on
    
    ts=Event(Clust==MSUnit(i))./sampling_rate; %in s
    move=ts(Moving(1)<ts & ts<Moving(2));
    static=ts(Static(1)<ts & ts<Static(2));
  
    subplot(2,3,1);
    %plot 1 firing rate over time
     binrange_static=[Static(1):1:Static(2)];
     histogram(static,binrange_static);
     title(strcat('Static FR of Unit ',num2str(MSUnit(i))));
     subplot(2,3,2);
    binrange_moving=[Moving(1):1:Moving(2)];
      histogram(move,binrange_moving);
       title(strcat('Moving FR of Unit ',num2str(MSUnit(i))));
      subplot(2,3,3) %moving vs staic line plot
      [static_count static_edges]=histcounts(static);
      [move_count move_edges]=histcounts(move);
      plot(static_edges(2:end),static_count,'LineWidth',3.5,'Color','b');
       hold on
        plot(move_edges(2:end),move_count,'LineWidth',3.5,'Color','r');
        title('FR of 2 states')
        legend('static','moving');
        hold off
      %box plot of distribution
      subplot(2,3,4)
      g=[ones(length(static_count'),1);2*ones(length(move_count'),1)];
       boxplot([static_count';move_count'],g,'Labels',{'Static','Moving'},'BoxStyle','outline','Colors','g');
       [h p]=ttest2(static_count',move_count');
       if h==0
         annotation('textbox', [0, 0.5, 0, 0], 'string', strcat('Static Moving state from same population ',' p=',num2str(p),' spontaneous'));
         neuron_type_stim(i,2)=1;
       elseif h==1 && mean(static_count)<mean(move_count)
         annotation('textbox', [0, 0.5, 0, 0], 'string', strcat('Static Moving state from different population ',' p=',num2str(p),' excitatory'));
         neuron_type_stim(i,2)=2;
       elseif h==1 && mean(static_count)>mean(move_count)
           annotation('textbox', [0, 0.5, 0, 0],'string', strcat('Static Moving state from different population ',' p=',num2str(p),' inhibitory'));
         neuron_type_stim(i,2)=3;
       end
       %ACG of two state
       [ISIm,ACGm] = ISI_ACG(move);
       [ISIs,ACGs] = ISI_ACG(static);
       range=[-50:0.5:50];
      subplot(2,3,5)%moving histogram
      histogram(ACGm,range);
      title('moving state ACG');
      xlabel('ACG(ms)');
      ylabel('Count');
      subplot(2,3,6)%static histogram
      histogram(ACGs,range);
      title('static state ACG');
      xlabel('ACG(ms)');
      ylabel('Count');
      saveas(f,fullfile(save_path,strcat('behavior state of unit ',num2str(MSUnit(i)),'.png')));
end
%COUNT matrix of stim neuron type
stim_type_count=zeros(3,1);
stim_type_ratio=zeros(3,1);
for c=1:3
stim_type_count(c,1)=sum(neuron_type_stim(:,2)==c); %1st row spon, 2nd row exci, 3rd row inhibitory
stim_type_ratio(c,1)=sum(neuron_type_stim(:,2)==c)/numel(MSUnit);
end








