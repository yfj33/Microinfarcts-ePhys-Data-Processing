function [ts_curation,label_old,label_new,firings_manual_curation] = MergedPair()
%Merged Pair Visualization: Show ISI, ACG merged 
%Input: read curation file
%Output: create a folder to store following  original cluster labels and
%new clusters label ,new firings.mat, all saved in the curation file


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

Clust=Firings(3,:);%cluster number
Event=Firings(2,:);%firing event occurence
PriCh=Firings(1,:);%primary channel
%Preallocate cells to store event time stamps for each cluster
Fir=cell(ClustN,1);
for c=1:ClustN
    Fir{c,1}=Event(Clust==ClustList(c));
end

pairs=[];
 reduction=0; %count how many clusters should be got rid of
 for n=1:numel(Merged_Pairs)

     reduction=reduction+[numel(str2num(Merged_Pairs{n}))-1];
    pairs=[pairs str2num(Merged_Pairs{n})];
 end
cluster_remain=MSUnit(~ismember(MSUnit,pairs));
new_dimen=numel(MSUnit)-reduction;

%preallocate new firings.mda
Fir_new=cell(numel(MSUnit),1);
%preallocate old cluster label
label_old=zeros(numel(MSUnit),1);
label_new=zeros(new_dimen,1);

for n=1:numel(Merged_Pairs)
  candi=str2num(Merged_Pairs{n});
  clus=[];
 
     for m=1:numel(candi) %each cell
         cluslabel=candi(m);
         clus=[clus Fir{cluslabel,1}];
  
     end

     Fir_new{find(MSUnit==min(candi)),1}=clus; %put merged cluster at the first cluster to be merged
     label_old(min(candi),1)=min(candi);
end
%fill in remaining unchanged cluster and their label
for m=1:numel(cluster_remain)
    Fir_new{find(cluster_remain(m)==MSUnit),1}=Fir{cluster_remain(m),1};
    label_old(cluster_remain(m),1)=cluster_remain(m);
end
%get rid of empty cell
ts_curation=Fir_new(~cellfun(@isempty,Fir_new));
save(fullfile(pwd,foldername,'ts_curation.mat'),'ts_curation','-mat');
%old label with merged clusters first cluster as label number
label_old(label_old==0)=[]; %old arrangement of label number
label_new=1:1:numel(ts_curation); %new arrangement of label number
save(fullfile(pwd,foldername,'label_new.mat'),"label_new",'-mat');
save(fullfile(pwd,foldername,'label_old.mat'),"label_old",'-mat');

%write friings.m file after merging, adjust at the original firing cluster

%delete rejected units
Goaway_unit=[Noise;Positive_Spikes];
event_manualcure=Event;
clust_manualcure=Clust;
for p=1:numel(Goaway_unit)
event_manualcure(Clust==Goaway_unit(p))=nan;

clust_manualcure(Clust==Goaway_unit(p))=nan;
end
event_manualcure(isnan(event_manualcure))=[];
clust_manualcure(isnan(clust_manualcure))=[];

%combine units by assigning new label to combined units
for k=1:numel(Merged_Pairs)
    candi=str2num(Merged_Pairs{k});
    assigned_lb=min(candi);
    for l=1:numel(candi)
    clust_manualcure(clust_manualcure==candi(l))=assigned_lb;
    end
end

firings_manual_curation=[event_manualcure;clust_manualcure];
save(fullfile(pwd,foldername,'firings_manual_curation.mat'),"firings_manual_curation",'-mat');

end
    






