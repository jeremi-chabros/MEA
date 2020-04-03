% script that takes a fraction of a voltage trace and plots before and
% after ttx

clear all; close all
filename                    =   '191210_FTD_slice1_DIV_g04_2018.mat';
prop_to_cut                 =   0.5; %proportion of recording to take
load(filename)

dat=dat(1:end*prop_to_cut,:);
%overwrite file
disp('saving data...')
save(strcat(filename(1:end-4),'_edited.mat'), 'ADCz','dat','channels','fs','header','uV','-v7.3');

clear all
