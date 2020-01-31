%script to run through and get filtered data from all files

%based on detectspikes.m Author: Tim Sit, sitpakhang@gmail.com 

clear all

%% load data
    %input data needs to be n samples x n channels
    %cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\mats'
    files=dir('*190705*.mat');
    files = files(~contains({files.name}, 'Spikes'));%remove unwanted files
    progressbar

        %% Filter signal 
        for file=1:length(files)
                
                %cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\filt_mats'
                cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
                if ~exist(strcat(files(file).name(1:end-4), '_Filtd', '.mat'))... %if file already exists, move on
                & files(file).name(end-8) ~= 'F' %if filename contain _Filtd this will be 0 and so will move on
                
                %cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\mats'
                               
                tic;
                load(files(file).name,'dat','channels','fs');
                disp(strcat('loaded: ',files(file).name))
  
                %filter         
                lowpass = 600; 
                highpass = 8000; 
                wn = [lowpass highpass] / (fs / 2); 
                filterOrder = 3; %changed from 3 to 5 by AD; and back
                [b, a] = butter(filterOrder, wn); 
                filteredMatrix = filtfilt(b, a, double(dat)); 
    
                %save as .mat matrix n sample x n channel
                fileName = strcat(files(file).name(1:end-4), '_Filtd', '.mat'); 
                % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
                disp('saving...')
                save(fileName, 'filteredMatrix','channels','fs', '-v7.3');
                    
                toc
                
                
                clear dat
                clear filteredMatrix
                clear channels
               else
                   disp(strcat(files(file).name(1:end-4),'_file done already'))
               end
               progressbar(file/length(files));
        end
    
