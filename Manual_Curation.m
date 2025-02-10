%% Sort figures into categories: all clusters noise, clusters(sub category narrow, wide, pyramidal)
clc;clear;close all
animal_folder=pwd;
folders=dir;
folders(1:2)=[];
folders_name={folders.name};
desired_folder=cellfun(@(y) y(1)=='x',folders_name);
desired_list=folders_name(desired_folder);


for i=1:numel(desired_list)

cd(fullfile(pwd,desired_list{i}));
templates=readmda('templates.mda');
cell_type=load('cell_type');
Combine=cell_type.cell_type.Combine; %1*N row for All Units
unit_cell=Combine(:,2);
MSUnit=Combine(:,1);
unit_num=size(templates,3);
AllClust=[1:1:unit_num];
% noise=xlsread(fullfile(session_name{n},'curation.xlsx'),'A:A');
% positive_spikes=xlsread(fullfile(session_name{n},'curation.xlsx'),'B:B');
% Noise=rmmissing(unique(noise));
% Positive_Spikes=rmmissing(unique(positive_spikes));
% Wasted=[Noise;Positive_Spikes]; %All useless units
% MSUnit=AllClust(~ismember(AllClust,[Noise;Positive_Spikes])); %All good units including Multi and Single units 
%Calculate FR before going to each day's folder
%for n=1:numel(Positive_Spikes)
%read cluster metrics 
%fid_met=fopen('cluster_metrics.json');
%met_info=fread(fid_met,inf);
%met_str=char(met_info');
%fclose(fid_met);
%Metrics=jsondecode(met_str);
%FR=FR+Metrics.clusters(Positive_Spikes(n)).metrics.firing_rate;

%move all units to All Cluster Folder
cufolder=fullfile(pwd,'clusters figures');%current folder where all figures is at
all_clusters=fullfile(cufolder,'All Clusters');
mkdir(all_clusters);
parfor a=1:unit_num
    cname=strcat('Cluster',num2str(AllClust(a)));
    
    movefile(fullfile(cufolder,strcat(cname,'.fig')),all_clusters);

    movefile(fullfile(cufolder,strcat(cname,'.png')),all_clusters);

end
%move units to folder Good Clusters folder with Subfolder for three cell
%types
good_clusters=fullfile(cufolder,'Good Clusters');
mkdir(good_clusters);
%make subfolder for three types of cell
narrow_folder=fullfile(good_clusters,'Narrow Interneurons');
wide_folder=fullfile(good_clusters,'Wide Interneurons');
pyramidal_folder=fullfile(good_clusters,'Pyramidal Cells');
mkdir(wide_folder);
mkdir(pyramidal_folder);
mkdir(narrow_folder);

for g=1:numel(MSUnit)
    fname=strcat('Cluster',num2str(MSUnit(g)),'.png');
    sourcefile=fullfile(all_clusters,fname);
    if unit_cell(g)==1;
    dest=narrow_folder;
    elseif unit_cell(g)==2;
      dest=wide_folder;
    elseif unit_cell(g)==3;
      dest=pyramidal_folder;
    end
  copyfile(sourcefile,dest);
end

%Read curation results from each column and masked out rejected units





% %go to figure folders
% cd(fullfile(pwd,'clusters figures'))
% clusters=dir('*.fig');
% %clusters(1:2)=[];
% %clusters_name=clusters.name;
% num=numel(clusters); %number of clusters
% AllClust=[1:1:num];
% AllClust=AllClust';
% %temp=AllClust(~ismember(AllClust,Noise));
% MSUnit=AllClust(~ismember(AllClust,[Noise;Positive_Spikes]));%All Units filtered out
% %put noise clusters in a folder
% noise_folder=fullfile(pwd,'noise clusters');
% mkdir(noise_folder);
% parfor i=1:numel(Noise)
% noise_cluster_name=strcat('Cluster',num2str(Noise(i)));
% %noise_fig=;
% %noise_png=strcat(noise_cluster_name,'.png');
% movefile(strcat(noise_cluster_name,'.fig'),noise_folder);
% movefile(strcat(noise_cluster_name,'.png'),noise_folder);
% %cd ..
% end
% 
% %put positive spikes in a folder
% pospike_folder=fullfile(pwd,'positive clusters');
% mkdir(pospike_folder);
% parfor i=1:numel(Positive_Spikes)
%     pospike_cluster_name=strcat('Cluster',num2str(Positive_Spikes(i)));
%    % pospike_fig=strcat(pospike,'.fig');
%    % pospike_png=strcat(pospike,'.png');
%     movefile(strcat(pospike_cluster_name,'.fig'),pospike_folder);
%     %cd ..
%     movefile(strcat(pospike_cluster_name,'.png'),pospike_folder);
%     %cd ..
% end
% 
% %the rest clusters is either a single/multi units
% Units_folder=fullfile(pwd,'negative units');
% mkdir(Units_folder);
% parfor i=1:numel(MSUnit)
%     ms_unit_name=strcat('Cluster',num2str(MSUnit(i)));
%    % ms_unit_fig=strcat(ms_unit_name,'.fig');
%    % ms_unit_png=strcat(ms_unit_name,'.png');
%     movefile(strcat(ms_unit_name,'.fig'),Units_folder);
%    % cd ..
%     movefile(strcat(ms_unit_name,'.png'),Units_folder);
%    % cd ..
% end
% 
% %negative cluster folder separare multiunit and single units
cd (animal_folder);

end
%save('Firing Rate','Firing_Rate','-mat');
%plot FR over days
%f=figure;
%plot([1:1:numel(week_info)],firing_rate,'FontSize','18');
%xticklabels=week_info;
%xlabel('Days');
%ylabel('Firing Rate (Hz)');
%saveas(f,'Firing Rate Plot');


