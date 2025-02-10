%compute the rejected cluster number and cell type ratio using Jiaao's npz
%file after adjusting parameters using old pipelines

%%read cluster rejection mask from folder of interest
folders=dir(pwd);
folders_name={folders.name};
FOI=cellfun(@(y) y(1)=='2',folders_name);
FOI_name=folders_name(FOI);

for n=1:numel(FOI_name)
animal_folder=fullfile(pwd,FOI_name{n});
mask_path=fullfile(animal_folder,'cluster_rejection_mask.mat');
rejection_mask=load(mask_path);
SU_mask=rejection_mask.single_unit_mask;
%CLUSTER NOISE SU
cluster_num=[1:1:numel(SU_mask)];
SU_cluster_num=cluster_num(SU_mask);
Noise_cluster_num=cluster_num(~ismember(cluster_num,SU_cluster_num));
%output xls file to the folder 
xlswrite(fullfile(animal_folder,'curation.xlsx'),Noise_cluster_num');
end
run MS_ChannelsView_V2.m
run cell_type_2022.m
run Move_Clusters2Folder.m











