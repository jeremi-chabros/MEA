    % script to cut down .mat data files to 6 mins ready for spike detection

    clear all
    data_and_scripts_dir = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output';
    cd(data_and_scripts_dir)
    files = dir('*FTD*.mat');                   % where your .mat files are
    files = files(~contains({files.name}, 'TTX'));
    files = files(~contains({files.name}, 'ttx'));
    files = files(~contains({files.name}, 'Spikes'));
    files = files(~contains({files.name}, 'stim'));
    files = files(~contains({files.name}, 'edited'));
    files = files(~contains({files.name}, '2min'));
    files = files(~contains({files.name}, '2001'));
    files = files(~contains({files.name}, '191210'));

    progressbar
    
    length_desired = 360; %time in secs 
    for i = 1:length(files)
            load(files(i).name)
        dat = dat(1:length_desired*fs,:);
        save(files(i).name, 'ADCz','dat','channels','fs','header','uV','-v7.3');
        progressbar(i/length(files))
    end