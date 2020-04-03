% combine figures (PNGs) produced by main script to a pdf
% create a folder and save the pdf and PNGs in there for each recording
%   may require moving the PNGs already created

% this script requires installing ghostscript:
% https://www.ghostscript.com/download.html
% this script also require the export_fig package
% https://uk.mathworks.com/matlabcentral/fileexchange/23629-export_fig

% inputs: files; spike method and parameter
clear all
files = dir('*FTD*ttx*.mat');
% files = files(~~contains({files.name}, 'TTX'));
files = files(~contains({files.name}, 'TTX2'));
files = files(~contains({files.name}, 'ttx4'));
files = files(~contains({files.name}, 'Spikes'));
files = files(~contains({files.name}, 'stim'));
files = files(~contains({files.name}, 'edited'));
files = files(~contains({files.name}, '2min'));
% files = files(~contains({files.name}, '2001'));
files = files(~contains({files.name}, '191210'));

method      = 'aSpikes'; %use mSpikes_, cSpikes_ or aSpikes (no _ if aSpikes)
parameter   = ''; % put value in inv. commas ('')

progressbar('files','figures')
for file = 1:length(files)
    filename = files(file).name;
    
    fig_names = dir(strcat(filename(1:end-4),'*.PNG'));
    fig_names = fig_names(~~contains({fig_names.name},strcat(method,parameter,'_')));
    fig_names = fig_names(~contains({fig_names.name},'based_on'));
    
    for fig = 1:length(fig_names)
        figure
        figname = fig_names(fig).name;
        imshow(imread(figname));
        
        h=gcf;
        set(h,'PaperPositionMode','auto');
        set(h,'PaperOrientation','landscape');
        set(h,'Position',[50 50 1200 800]);
        print(gcf, '-dpdf', strcat(figname(1:end-4),'.pdf'))
        
        close(h)
        
        progressbar(file/length(files),fig/length(fig_names))
    end
    
    %put pdf figures in list of cells
    fig_PDFnames = dir(strcat(filename(1:end-4),'*.pdf'));
    
    for figname = 1:length(fig_PDFnames)
        fig_PDFnames_list{figname,1} = fig_PDFnames(figname).name;
    end
    
    % append figures and save
    try
    append_pdfs(strcat(filename(1:end-4),'_',method,parameter,'_all_figs.pdf'),fig_PDFnames_list{:});
    catch
        disp('already done')
    end
    
    % create folder for the recording and move all files into that
    if ~exist(filename(1:end-4), 'dir')
        mkdir(filename(1:end-4))
    end
    
    files2move = dir(strcat(filename(1:end-4),'*.p*'));
    files2move = files2move(~~contains({files2move.name},strcat(method,parameter,'_')));
    
    cd(filename(1:end-4))
    if ~exist(strcat(filename(1:end-4),'_',method,parameter), 'dir')
        mkdir(strcat(filename(1:end-4),'_',method,parameter))
    end
    
    cd ..
    
    for file_m = 1:length(files2move)
        movefile(files2move(file_m).name,strcat(filename(1:end-4),'/',filename(1:end-4),'_',method,parameter));
    end
    

    progressbar(file/length(files),fig/length(fig_names))
        clear fig_names
end