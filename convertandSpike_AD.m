% convert files and get mSpikes; then c spikes
cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.1.File_Conversion_Scripts'
MEAbatchConvert_alex
cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
batchGetSpike %MS method
batchGetSpikeC %cwt method

