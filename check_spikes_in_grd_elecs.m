% last edit AD 16-10-2019

% this script should be run after grd_elecs.m to check that the grounded
% electrodes contain ~0 spikes.
% you need to specify the names of the files to check (which usually 
% depends on the spike detection parameters used)
% you also need to specify the name of the .csv file containing the names
% of the .mat files and the electrodes to ground — this should have been
% determined before running grd_elecs.m so should be the same here.

Es_to_grd_table=readcell('to_ground.csv');      %name of the .csv file 
%containing the table of names of the recordings and the electrodes to 
%ground in each of those recordings

files = dir('*190705*Spikes*L-0.1254*.mat');  % desired files to check
%files = files(~contains({files.name}, 'Spikes'));%remove unwanted files


progressbar
for i = 1:size(Es_to_grd_table,2)
    tic;
    
for j=2:length(Es_to_grd_table)
    Es_to_grd(j-1)=Es_to_grd_table{j,i};
end

load(files(i).name)

a=sum(full(cSpikes));
for k = 1:length(Es_to_grd)
    b(i,k)=a(find(channels==Es_to_grd(k)));
end


    toc
    progressbar(i/size(Es_to_grd_table,2));
end

b %each row is a different file (in the order listed in files variabel;
    %each column is a channel (in the order listed in Es_to_grd variable)