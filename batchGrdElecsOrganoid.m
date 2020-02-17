%200127_FTDOrg_GrpA_3A_Slice1_mSpikes_3_adjM.mat

%create cell array containing file names

clear all; close all
files1 = dir('*FTD*.mat');
files1 = files1(~~contains({files1.name}, 'mSpikes_3'));
files1 = files1(~contains({files1.name}, 'edited'));
files1 = files1(~contains({files1.name}, '2min'));
files1 = files1(~contains({files1.name}, 'adjM'));
files1 = files1(~contains({files1.name}, '191210'));
for i = 1:length(files1)
    files{i,1} = files1(i).name;
end

% now open the cell array and manually enter the IDs of the electrodes to
% ground
for i = 1:length(files1)
    files{i,2} = [85];
end

%190814
files{1,2} = [85,86];
files{2,2} = [85,86];
files{3,2} = [85,86];
files{4,2} = [85,86];
files{5,2} = [85,86];

%191209
files{6,2} = [85,86,17];
files{7,2} = [85,86];
files{8,2} = [85,86];
files{9,2} = [85,86];
files{10,2} = [85,86];
files{11,2} = [85,86];
files{12,2} = [85,86];
files{13,2} = [85,86];
files{14,2} = [85,86];
files{15,2} = [85,86];

%200114
files{16,2} = [85,86];
files{17,2} = [85,86];
files{18,2} = [85,86];
files{19,2} = [85,86,23];
files{20,2} = [85,86];
files{21,2} = [85,86];
files{22,2} = [85,86];
files{23,2} = [85,86];
files{24,2} = [85,86];
files{25,2} = [85,86];
files{26,2} = [85,86,56,55];
files{27,2} = [85,86];
files{28,2} = [85,86];
files{29,2} = [85,86];
files{30,2} = [85,86];
files{31,2} = [85,86];
files{32,2} = [85,86,56];

%200127
files{33,2} = [85,86,25];
files{34,2} = [85,86,23,31,54];
files{35,2} = [85,86];
files{36,2} = [85,86,53,54];
files{37,2} = [85,86,17];
files{38,2} = [85,86];
files{39,2} = [85,86,17];
files{40,2} = [17,85];


progressbar
for i = 1:length(files)
    filename = files{i,1};
    channelIDs = files{i,2};
    grd_elec_organoid_fcn(filename,channelIDs)
    progressbar(i/length(files))
end
