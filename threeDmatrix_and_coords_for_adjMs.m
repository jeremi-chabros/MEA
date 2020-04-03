%% load files and create 1 3d matrix

clear all

files = dir('*het*adjM*DA60*');
for j=1:length(files)
    load(files(j).name,'adjM');
    adjM_all(:,:,j)=adjM; 
    clear adjM
end

 
%% load channels and create 60 x 2 vector

 
    load(files(1).name);
    coordstr=int2str(channels);
    for i=1:length(coordstr)
        
    coord(i,1)=str2double(coordstr(i));
    coord(i,2)=str2double(coordstr(i+60));
    end
    
%% save
fileName=strcat('3d_matrix_coords_MEAs_het')
save(fileName,'adjM_all','coord')

