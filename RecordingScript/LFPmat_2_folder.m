function [subfolderlist] = LFPmat_2_folder(current_folder)
%move LFP.mat file from each session to a folder called LFPs
cd(current_folder);
mkdir('Stroke LFPs');
subfolderlist=dir;
subfoldername={subfolderlist.name};
desiredfolder = cellfun(@(y) y(1)=='2' , subfoldername );
folderlist=subfoldername(desiredfolder); % foldername of each session

for i=1:numel(folderlist)
    cd(folderlist{i})
    matinfo=dir('*.mat')
    matname={matinfo.name}
    movefile(char(matname),[current_folder '\' 'Stroke LFPs']);
    cd(current_folder)
    
end

end

