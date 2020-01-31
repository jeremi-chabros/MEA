% Last update: AD 15-10-2019
%GROUNDING ELECTRODES FUNCTION
%{ 
inputs:     1. row vector of electrode IDs that you want to ground
            this should be in the MCS format (column number
            followed by row number ##)

            2. variable called channels in which the channel IDs 
            in MCS format are list in the order that they appear
            in the data file (dat variable) containing the voltage
            traces.
        
            3. variable containing the voltage traces
            e.g. for a 12 min recording at 25 kHz, this will be an
            18,000,000 * 60 matrix
            
            4. filename that you want to overwrite (including the
            file extension .mat)

outputs:    This will output the new data file with the chosen electrodes
            grounded (set to equate to the reference electrode) and the 
            old data file will be overwritten.
%}

%function dat = grd_elecs(Es_to_grd,channel_list,data,filename)
%dat(:,find(channel_list==Es_to_grd))=dat(:,find(channel_list==15));
%end

%% as a stand alone script rather than a function...

%load table of filenames and electrodes to ground
%needs to be: top row = filenames; 
%                each subsequent row is an electrode MCS ID
Es_to_grd_table=readcell('to_ground.csv');

progressbar
for i = 1:size(Es_to_grd_table,2)
    tic;
%extract filename
current_file = Es_to_grd_table{1,i}

%create electrode ID vector 
for j=2:length(Es_to_grd_table)
    Es_to_grd(j-1)=Es_to_grd_table{j,i};
end
    
%ground electrodes:
load(current_file);
for h=1:length(Es_to_grd)
dat(:,find(channels==Es_to_grd(h)))=dat(:,find(channels==15));
end

%overwrite file
disp('saving data...')
save(current_file, 'ADCz','dat','channels','fs','header','uV','-v7.3');
clear dat channels ADCz ans current_file Es_to_grd fs header uV
    toc
    progressbar(i/size(Es_to_grd_table,2));
end

clear all




