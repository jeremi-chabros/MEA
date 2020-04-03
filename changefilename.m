%% script to rename tristan files to MPTYYMMDD_ID_DIV##
clear all

cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'

files=dir('191210_FTD*');

remove = 'FTD_';

progressbar('filenames changed')
for file=1:length(files)
    
    rm_index = [strfind(files(file).name,remove) : ...
        strfind(files(file).name,remove) - 1 + length(remove)];
    
    newname = files(file).name;
    newname(rm_index) = [];
    
    movefile(files(file).name,newname);
    progressbar(file/length(files))
end

%% mistake add two _s

%files=dir('*__*')
%movefile(files(1).name,strcat(files(1).name(1:13),files(1).name(15:end)));