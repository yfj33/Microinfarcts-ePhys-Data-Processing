function [spikes,spikes_single,results_merging,results_merging_SU] = MergedPair_V2(SR)
%[spikes,spikes_single,results_merging,results_merging_SU] = MergedPair_V2(SR)
%Input: SR
%Output: firings_manual_curation.mat, template_manual_curation.mat after
%merging, old and new labels after merging for MS units and single units;
%also spikes.mat for cellexplorer input


foldername='manual curation results'; %folder to save manual curation results
if ~exist(foldername)
    mkdir(foldername);
end
%Read curation results from each column and masked out rejected units
templates=readmda('templates.mda');
unit_num=size(templates,3);
AllClust=[1:1:unit_num];
noise=xlsread('curation.xlsx','A:A');
positive_spikes=xlsread('curation.xlsx','B:B');
Noise=rmmissing(unique(noise));
Positive_Spikes=rmmissing(unique(positive_spikes));
MSUnit=AllClust(~ismember(AllClust,[Noise;Positive_Spikes]));
%read curation pairs from excel file
Curation=readtable('curation.xlsx');
merged_pairs=Curation{:,3};
%Merged_pairs=str2num(merged_pairs);
Merged_Pairs=merged_pairs(~cellfun(@isempty,merged_pairs)); %in a cell with each merged pairs in a cluster
%read firings.mda and assign events to each cell
Firings=readmda('firings.mda');
ClustN=numel(unique(Firings(3,:)));%cluster number
ClustList=unique(Firings(3,:));%Cluster List based on sequence small to large

Clust=Firings(3,:);%cluster number/label
Event=Firings(2,:);%firing event occurence
PriCh=Firings(1,:);%primary channel
pairs=[];
 reduction=0; %count how many clusters should be got rid of
 for n=1:numel(Merged_Pairs)

     reduction=reduction+[numel(str2num(Merged_Pairs{n}))-1];
    pairs=[pairs str2num(Merged_Pairs{n})];
 end
cluster_remain=MSUnit(~ismember(MSUnit,pairs));
new_dimen=numel(MSUnit)-reduction;


%delete rejected units
Goaway_unit=[Noise;Positive_Spikes];
event_manualcure=Event;
clust_manualcure=Clust;
templates_merging=templates;
for p=1:numel(Goaway_unit)
event_manualcure(Clust==Goaway_unit(p))=nan;

clust_manualcure(Clust==Goaway_unit(p))=nan;

templates_merging(:,:,Goaway_unit(p))=nan;
end
event_manualcure(isnan(event_manualcure))=[];
clust_manualcure(isnan(clust_manualcure))=[];


%combine units by assigning new label to combined units: smallest label
%number remains and others deleted
for k=1:numel(Merged_Pairs)
    candi=str2num(Merged_Pairs{k});
    assigned_lb=min(candi);
    for l=1:numel(candi)
   
    clust_manualcure(clust_manualcure==candi(l))=assigned_lb;
    

    end

end
%new filtwaveform using current template by merging the template waveform
New_ClustList=unique(clust_manualcure);
New_Clustnum=numel(New_ClustList);
label_native=New_ClustList;%native neuron labels from MS 
label_new=1:1:New_Clustnum;

%calculate events number before merging labels
%times in secconds
times_b4=cell(1,numel(MSUnit));
Event_sec=Event/SR;
for i=1:numel(MSUnit)
times_b4{1,i}=Event_sec(Clust==MSUnit(i));  %1*N cell structure for N units each has M spikes
end
%total
total_b4=zeros(1,numel(MSUnit));
for i=1:numel(MSUnit)
    total_b4(1,i)=numel(times_b4{i});
end

%template after curation, preallocate
clip_num=100;
Templates_Merging=zeros(size(templates,1),clip_num,New_Clustnum);
for a=1:numel(Merged_Pairs)
     candi=str2num(Merged_Pairs{a});
     lb=min(candi);
     lb_new=label_new(lb==label_native);
     %wt=zeros(1,candi)%store number of events each time f
    
   mergepair=templates_merging(:,:,candi);
    wt=zeros(1,numel(candi));%events number
   for b=1:numel(candi)
      
       wt(b)=total_b4(candi(b)==MSUnit);%sequence messed up
       
   end
   temp=zeros(size(mergepair,1),size(mergepair,2),1);%
   for c=1:numel(candi)
       temp=mergepair(:,:,c).*wt(c)+temp;
   end
  Templates_Merging(:,:,lb_new)=temp/sum(wt);
   
end

%templates_merging(isnan(templates_merging))=[];% delete rejected clusters 
%Fill in the rest of the unchanged ones
for m=1:numel(cluster_remain)
    seq=label_new(cluster_remain(m)==label_native);
    Templates_Merging(:,:,seq)=templates_merging(:,:,cluster_remain(m));
end

%Detect Primary channel
PC_curation=zeros(1,New_Clustnum);
for i=1:New_Clustnum
    valley=min(Templates_Merging(:,:,i),[],2); %lowest point n*1 column
    [M id]=min(valley); %id is the channel number absolute, need to refer back to channel order
    PC_curation(1,i)=id; %PC of each cluster
    
end
%filtWaveform
filtWaveform=cell(1,New_Clustnum);
timeWaveform=cell(1,New_Clustnum);
for i=1:New_Clustnum
    %pc=unique(PC(MSUnit(i)==Clust));
    %filtWaveform{i}=templates(pc,:,MSUnit(i));
    filtWaveform{i}=Templates_Merging(PC_curation(i),:,i);
    timeWaveform{i}=[-49:1:50]/SR*1000; %unit in ms
end

% firings_manual_curation=[event_manualcure;clust_manualcure];
%save(fullfile(pwd,foldername,'firings_manual_curation.mat'),"firings_manual_curation",'-mat');


%primary channel
PC=zeros(1,numel(event_manualcure));
for z=1:numel(event_manualcure)
    PC(1,z)=PC_curation(clust_manualcure(z)==label_native);
end
%save(fullfile(pwd,foldername,'firings_after_merging.mat'),'Firings_Merging','-mat');
%save(fullfile(pwd,foldername,'label_native.mat'),'label_native','-mat');
%save(fullfile(pwd,foldername,'label_new.mat'),'label_new','-mat');
%save(fullfile(pwd,foldername,'templates_after_merging.mat'),'Templates_Merging','-mat');

%spikes 
%ts before curation and event number for each cluster
ts=cell(1,New_Clustnum);
times=cell(1,New_Clustnum);
total=zeros(1,New_Clustnum);
cluID=label_native;
shankID=ones(1,New_Clustnum);
spindices=[event_manualcure;clust_manualcure];
numcells=New_Clustnum;

for x=1:New_Clustnum
    ts{1,x}=event_manualcure(clust_manualcure==New_ClustList(x));
    times{1,x}= ts{1,x}./SR;
    total(1,x)=numel(ts{1,x});
end


%single unit mask
% single_unit_mask=ones(1,New_Clustnum);
% for i=1:numel(ts)
%     isi=diff(times{i});
%     if sum(isi<0.002)/numel(isi)>0.01
%         single_unit_mask(1,i)=0;
%     end
% 
% end
% 
% ts(find(single_unit_mask==0))=[];
% times(find(single_unit_mask==0))=[];
% cluID(find(single_unit_mask==0))=[];
% shankID(find(single_unit_mask==0))=[];
% filtWaveform(:,:,find(single_unit_mask==0))=[];
% timeWaveform(:,:,find(single_unit_mask==0))=[];
% lb_off=label_native;
% lb_off(1,find(single_unit_mask==0))=[];
% 
% event_manualcure(clust_manualcure==lb_off)=[];
% clust_manualcure(clust_manualcure==lb_off)=[];
% spindices=[event_manualcure;clust_manualcure];
% numcells=numel(times);
% for x=1:numel(ts)
%    
%     total(1,x)=numel(ts{1,x});
% end

%Compute only single units
spikes.ts=ts;
spikes.times=times;%in s
spikes.cluID=cluID;
spikes.shankID=shankID;
spikes.filtWaveform=filtWaveform;
spikes.timeWaveform=timeWaveform;
spikes.spindices=spindices;
spikes.numcells=numcells;
spikes.total=total;
%Firings after merging
Firings_Merging=[PC;event_manualcure;clust_manualcure];
results_merging.Firings_Merging=Firings_Merging;
results_merging.label_native=label_native;
results_merging.label_new=label_new;
results_merging.Templates_Merging=Templates_Merging;
results_merging.spikes=spikes;
save(fullfile(pwd,foldername,'results_merging.mat'),'results_merging','-mat');

%%compute single units only with merging results
%Templates_Merging_SU=Templates_Merging;
Firings_Merging_SU=Firings_Merging;
times_SU=times;
%single units label list
SU_mask=ones(1,numel(label_native));
for s=1:numel(SU_mask)
    isi=diff(times_SU{s});
    if sum(isi<0.002)/numel(isi)>0.01 %2% violation?
        SU_mask(s)=0;
    end
end
label_native_SU=label_native;
label_native_SU=label_native_SU(SU_mask==1);%cutoff multiunits isi 2ms >1%
label_new_SU=1:1:numel(label_native_SU);
label_native_MU=label_native(~ismember(label_native,label_native_SU));

%update SU firings and templates
%firings
for i=1:numel(label_native_MU)
   Firings_Merging_SU(:,find(label_native_MU(i)==Firings_Merging_SU(3,:)))=[];%each loop deletes the given label
end
%templates
 Templates_Merging_SU=zeros(size(Templates_Merging,1),size(Templates_Merging,2),numel(label_new_SU));

for k=1:numel(label_native_SU)
    idx=label_new(label_native_SU(k)==label_native);
    Templates_Merging_SU(:,:,k)=Templates_Merging(:,:,idx);
end

ts_SU=cell(1,numel(label_native_SU));
total_SU=zeros(1,numel(label_native_SU));
cluID_SU=label_native_SU;
shankID_SU=ones(1,numel(label_native_SU));
spindices_SU=Firings_Merging_SU(2:3,:);
numcells_SU=numel(label_native_SU);

filtWaveform_SU=cell(1,numcells_SU);
timeWaveform_SU=cell(1,numcells_SU);
Times_SU=cell(1,numcells_SU);
for j=1:numcells_SU

   % ts_SU{1,j}=times_SU{1,j}*SR;
  %  total_SU(1,j)=numel(ts_SU{1,j});
    temp_id=label_new(label_native_SU(j)==label_native);
    filtWaveform_SU{j}=filtWaveform{temp_id};
    timeWaveform_SU{j}=timeWaveform{1,j};
    Times_SU{1,j}=times{1,temp_id};
    ts_SU{1,j}= Times_SU{1,j}*SR;
    total_SU(1,j)=numel(ts_SU);
end

spikes_single.ts=ts_SU;
spikes_single.times=Times_SU;
spikes_single.cluID=cluID_SU;
spikes_single.filtWaveform= filtWaveform_SU;
spikes_single.timeWaveform=timeWaveform_SU;
spikes_single.spindices=spindices_SU;
spikes_single.numcells=numcells_SU;
spikes_single.total=total_SU;
spikes_single.shankID=shankID_SU;

%save single units results
results_merging_SU.Firings_Merging_SU= Firings_Merging_SU;
results_merging_SU.label_native_SU=label_native_SU;
results_merging_SU.label_new_SU=label_new_SU;
results_merging_SU.Templates_Merging_SU=Templates_Merging_SU;
results_merging_SU.spikes_single=spikes_single;
save(fullfile(pwd,foldername,'results_merging_SU.mat'),'results_merging_SU','-mat');

end








