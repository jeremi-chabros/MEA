close all
clear all
datadir = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output';
cd(datadir);
files = dir('*MPT*DIV28*Spikes*adjM*DA60*');
files = files(~contains({files.name}, 'het'));%remove unwanted files
for file=1:length(files)
    load(files(file).name,'adjM');
    figure
    imagesc(adjM)
end