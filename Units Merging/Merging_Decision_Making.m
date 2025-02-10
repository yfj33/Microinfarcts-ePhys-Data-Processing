%CCG Merging
%Compute CCG of two units and decide if they should be merged
clc;clear;close all
binSize=0.001; %in s
duration=1; %in s
sr=30000;
Units_2_Decide=[43 42];
celltype=load('cell_type_2023.mat');
Comb_SU=celltype.cell_type.Comb_SU;
spikes_SU=celltype.cell_type.spikes_SU;
waveform=spikes_SU.filtWaveform;
%preallocate for CCG computations
Groups=[];
Times=[];
WF=cell(1,numel(Units_2_Decide));
FR_times=[];
%assemble into one group for CCG computation
for i=1:numel(Units_2_Decide)
num=Units_2_Decide(i);
seq=find(num==Comb_SU(:,1))
Times=[Times spikes_SU.times{seq}];
FR_times=[FR_times numel(spikes_SU.times{seq})];
Groups=[Groups i-1+ones(1,numel(spikes_SU.times{seq}))];
WF{1,i}=waveform{seq};
end
%compuite ISI 
ISI=diff(Times);
ISI_2pct=sum(ISI<0.002)/numel(ISI);
fprintf(strcat('2ms Violation value is ',' ',num2str(ISI_2pct)));
if ISI_2pct<0.02
    fprintf('Merged spike train is still single unit');
else
    fprintf('Merged spike train does not meet single unit criteria');
end
[ccg,t]=CCG(Times,Groups,'binSize',binSize,'duration',duration,'Fs',1/sr);

%merged ACG
%OUTPUT
%   ccg     [t x ngroups x ngroups] matrix where ccg(t,i,j) is the
%           number (or rate) of events of group j at time lag t with  
%           respect to reference events from group i
%   t       time lag vector (units: seonds)
% Burst index is determined by calculating the average number of spikes in the 3-5 ms bins of the spike
%    autocorrelogram divided by the average number of spikes in the 200-300 ms bins.
[acg,t2]=CCG(Times,ones(1,numel(Times)),'binSize',binSize,'duration',duration,'Fs',1/sr);
BurstIndex_Royer2012 = mean(acg(binSize+1+3:binSize+1+5,1))/mean(acg(binSize+1+200:binSize+1+300,1));

%output ccg ccg  [t x ngroups x ngroups] matrix where ccg(t,i,j) is the
%number (or rate) of events of group j at time lag t with  
%respect to reference events from group i
%time lag vector (units: seonds)
%calculate merge unit waveform
Times=[-49:50]/30000;
avg_WF=(WF{1,1}*FR_times(1)+WF{1,2}*FR_times(2))/(FR_times(1)+FR_times(2));
[t_value t_id]=min(avg_WF);
[p_value p_id]=max(avg_WF)
T2P_avg=abs(Times(t_id)-Times(p_id));
ampt2p=p_value-t_value;
ampabs=abs(t_value);

figure
subplot(3,2,1)
bar(t*1000,ccg(:,2,1),'FaceColor','b','EdgeColor','none');
xlabel('Time(ms)');
title('Unit 2 as Reference CCG');
subplot(3,2,2)
bar(t*1000,ccg(:,1,2),'FaceColor','b','EdgeColor','none');
title('Unit 1 as Reference CCG');
subplot(3,2,3)
bar(t*1000,ccg(:,1,1),'FaceColor','b','EdgeColor','none');
title('Unit 1 ACG');
subplot(3,2,4)
bar(t*1000,ccg(:,2,2),'FaceColor','b','EdgeColor','none');
title('Unit 2 ACG');
subplot(3,2,5)
histogram(ISI*1000,'BinEdges',[0:0.5:500],'FaceColor','b','EdgeColor','none','Normalization','probability');
hold on
xline(2,'LineWidth',1,'color','r')
xlabel('Time(ms)');
subplot(3,2,6)
hold on
axis off
plot([-49:50]/30000,WF{1,1},'Color','r','LineWidth',2);
plot([-49:50]/30000,WF{1,2},'Color','b','LineWidth',2);
plot([-49:50]/30000,avg_WF,'Color','c','LineWidth',2);
text(-20/30000,-20,strcat('T2P merge ',num2str(T2P_avg*1000),'ms'));
text(-20/30000,-40,strcat('Bursting Index ',num2str(BurstIndex_Royer2012)));
text(-20/30000,-60,strcat('ampt2p ',num2str(ampt2p)));
text(-20/30000,-80,strcat('ampabs ',num2str(ampabs)));
legend(num2str(Units_2_Decide(1)),num2str(Units_2_Decide(2)),'merging waveform');

% merging_info=
% save(fullfile(pwd,'Merging Info'),'merging_info');








