% loop to correct NaNs
for check=1:length(output)
    
    if isnan(output(check).netw_density)
        output(check).netw_density=0;
    else
    end

    
end