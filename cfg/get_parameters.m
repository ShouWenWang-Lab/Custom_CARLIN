function params = get_parameters(cfg)

    [folder, ~, ~] = fileparts(mfilename('fullpath'));
    
    params = inputParser;
    params.KeepUnmatched = false;
    params.CaseSensitive=true;
    params.PartialMatching=false;
    
    addRequired(params,'fastq_file');
    addRequired(params,'cfg_type');
    addRequired(params,'outdir');
    
    addParameter(params,'read_cutoff_UMI_denoised',10);
    addParameter(params,'read_override_UMI_denoised', NaN);
    
    if (strcmp(cfg.type, 'SC'))
        
        addParameter(params,'max_cells',inf);        
        addParameter(params,'read_cutoff_CB_denoised',10);
        addParameter(params,'read_override_CB_denoised', NaN);
        
        if (strcmp(cfg.SC.Platform, 'InDrops'))            
            if (cfg.SC.Version==2)
                addParameter(params, 'ref_CB_file', sprintf('%s/InDropsV2_barcodes.txt.gz', folder));
            elseif (cfg.SC.Version==3)
                addParameter(params, 'ref_CB_file', sprintf('%s/InDropsV3_barcodes.txt.gz', folder));
            else
                error('Unsupported InDrops version');
            end
        elseif (strcmp(cfg.SC.Platform, '10x'))
            if (cfg.SC.Version==2)
                addParameter(params, 'ref_CB_file', sprintf('%s/10xV2_barcodes.txt.gz', folder));            
            elseif (cfg.SC.Version==3)
                addParameter(params, 'ref_CB_file', sprintf('%s/10xV3_barcodes.txt.gz', folder));
            else
                error('Unsupported 10x version');
            end
        else
            error('Unsupported SingleCell platform');
        end
    else
        addParameter(params,'max_molecules',inf);
    end
    
end