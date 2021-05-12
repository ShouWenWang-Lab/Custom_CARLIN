function analyze_CARLIN(fastq_file, cfg_type, outdir, varargin)
%analyze_CARLIN calls alleles from FASTQs sequencing the CARLIN amplicon. 
%
%   analyze_CARLIN(FASTQ_FILE, CFG_TYPE, OUTDIR) analyzes FASTQ_FILE 
%   (*.fastq, *.fastq.gz) generated according to CFG_TYPE (one of 'BulkDNA', 
%   'BulkRNA', 'scInDropsV2', 'scInDropsV3', 'sc10xV2', 'sc10xV3') and saves the 
%   following files to OUTDIR:
%
%   Analysis.mat.gz - contains all variables including a compactified 
%     representation of the FastQ files, a depot of unique aligned
%     sequences, the collection of CBs and UMIs detected before and after 
%     denoising, and the called alleles. This file may be large, and is not
%     intended for use by a casual-user, rather to save state to do any
%     custom one-off analysis.
%
%   Summary.mat - contains a subset of the variables in Analysis.mat that 
%     are usually sufficient for subsequent downstream analysis in MATLAB
%     including a list of alleles, their frequencies, tags (UMIs or CBs) 
%     reporting that allele, input parameters, thresholds used by the 
%     algorithm - and, for SC runs - the reference barcode list.
%
%   Alleles.png - Plot of sequence and distribution of alleles with
%        summary statistics.
%
%   AlleleAnnotations.txt - Annotations describing each allele in plaintext
%     format (HGVS) for use by other downstream tools.
%
%   AlleleColonies.txt - List of tags (UMIs or CBs) which report each
%     allele, for use by other downstream tools. The ordering corresponds 
%     to AlleleAnnotations.txt
%
%   Results.txt - Human-readable summary of the analysis results.
%
%   Diagnostic.png - Detailed dashboard of various diagnostic quantities.
%
%   Warnings.txt - List of issues detected in the data.
%
%   Log.txt - Running log of the pipeline.
%
%   The pipeline can optionally be invoked with some extra paramters:
%
%   analyze_CARLIN(..., 'read_cutoff_UMI_denoised', cutoff) uses a minimum
%   read threshold of 'cutoff' when attempting to call alleles from denoised 
%   UMIs (default = 10). This cutoff is not used in its raw form, but
%   combined as part of a cutoff function. It represents an absolute floor,
%   but the cutoff used in practice will generally be higher.
%
%   analyze_CARLIN(..., 'read_override_UMI_denoised', cutoff) short circuits
%   the cutoff function, and sets the threshold to 'cutoff'. All UMIs with a
%   read count >= 'read_override_UMI_denoised' will be asked to call an allele. 
%   Default is unset.
%
%   For bulk sequencing runs (CFG_TYPE='Bulk*'):
%
%   analyze_CARLIN(..., 'max_molecules', N) considers (at most) the N denoised
%   UMIs with the most reads, for calling alleles. For CFG_TYPE='BulkDNA',
%   N is the number of cells in the sample. For CFG_TYPE='BulkRNA', N is the 
%   number of transcripts (if unsure, 10x the number of cells is suitable
%   guess). Default is inf.
%   
%   For single-cell sequencing runs of the CARLIN amplicon (CFG_TYPE='sc*'):
%
%   analyze_CARLIN(..., 'max_cells', N) considers (at most) the N denoised
%   cell barcodes with the most reads, for calling alleles. Defaults to 
%   the number of possible CBs in the single-cell platform specified by
%   CFG_TYPE.
%
%   analyze_CARLIN(..., 'read_cutoff_CB_denoised', cutoff) uses a minimum
%   read threshold of 'cutoff' when attempting to call alleles from denoised 
%   CBs (default = 10). This cutoff is not used in its raw form, but
%   combined as part of a cutoff function. It represents an absolute floor,
%   but the cutoff used in practice will generally be higher.
%
%   analyze_CARLIN(..., 'read_override_CB_denoised', cutoff) short circuits
%   the cutoff function, and sets the threshold to 'cutoff'. All CBs with a
%   read count >= 'read_override_CB_denoised' will be asked to call an allele. 
%   Default is unset.
%
%   analyze_CARLIN(..., 'ref_CB_file', file) uses the reference list of 
%   cell barcodes in the file specified by 'ref_CB_file' when denoising 
%   barcodes found in the FastQ files. The reference list should have one 
%   cell barcode per line. Each cell barcode should be a string consisting 
%   of only the characters {A,C,G,T}. The length of the barcode should match
%   the length expected by the platform specified by CFG_TYPE. This reference
%   list is typically produced by the software used to process the 
%   corresponding transcriptome run. Defaults to the full barcode list in 
%   the single-cell platform specified by CFG_TYPE. Specifying a reference
%   list leads to more accurate denoising than using the full barcode list.
%
%   Examples:
%
%       % Run analysis on bulk RNA sample
%       analyze_CARLIN('run1/PE.fastq.gz', 'BulkRNA', 'run1_output');
%
%       % Run analysis on multiple sequencing runs of the same bulk DNA sample
%       analyze_CARLIN({'run1/PE.fastq.gz'; run2/PE.fastq.gz'}, 'BulkDNA', 'combined_output');
%
%       % Run analysis on sample prepared with InDropsV3
%       analyze_CARLIN('indrops/amplicon/Lib1.fastq', 'scInDropsV3', ...
%                      'output', 'ref_CB_file', 'indrops/transcriptome/abundant_barcodes.txt');
%
%       % Run analysis on sample prepared with 10xGenomics V2
%       analyze_CARLIN({'tenx/amplicon/Mouse_R1_001.fastq', 'tenx/amplicon/Mouse_R2_001.fastq'}, ...
%                      'sc10xV2', 'tenx/processed_amplicon', ...
%                      'ref_CB_file', 'tenx/transcriptome/filtered_barcodes_umi_mt.txt');
% 
%   If you use this code, please cite:
%
%   S. Bowling, D. Sritharan, F. G. Osorio, M. Nguyen, P. Cheung, 
%   A. Rodiguez-Fraticelli, S. Patel, W-C. Yuan, Y. Fujiwara, B. E. Li, S. H. Orkin, 
%   S. Hormoz, F. D. Camargo. "An Engineered CRISPR-Cas9 Mouse Line for 
%   Simultaneous Readout of Lineage Histories and Gene Expression Profiles 
%   in Single Cells." Cell (2020), https://doi.org/10.1016/j.cell.2020.04.048
%
%   Author: Duluxan Sritharan. Hormoz Lab. Harvard Medical School.

    tic;

    diary([tempname '.txt']);
    diary on;
    fprintf('Writing log to temporary location: %s\n', get(0,'DiaryFile'));
    
    % Parse parameters
    cfg = parse_config_file(cfg_type);
    params = get_parameters(cfg);    
    parse(params, fastq_file, cfg_type, outdir, varargin{:});
    clear fastq_file cfg_type outdir varargin
        
    % Setup directory and log file
    if (~exist(params.Results.outdir, 'dir'))
        mkdir(params.Results.outdir);        
    end

    % Make FQ representation
    if (strcmp(cfg.type, 'Bulk'))    
        FQ = BulkFastQData(params.Results.fastq_file, cfg);
    else
        ref_CBs = get_SC_ref_BCs(params.Results.ref_CB_file);
        FQ = SCFastQData(params.Results.fastq_file, cfg);                         
    end

    % Align unique sequences
    aligned = AlignedSEQDepot(FQ.get_SEQs());
    aligned.sanitize_prefix_postfix();
    aligned.sanitize_conserved_regions();
    
    % Make CB collection and call alleles
    if (strcmp(cfg.type, 'Bulk'))
        tag_collection = BulkUMICollection.FromFQ(FQ);
        [tag_collection_denoised, tag_denoise_map] = tag_collection.denoise(aligned);
        %thresholds=struct('one_tenth_99_pctl', 1, 'max_molecules', 1, 'equal_partition', 0, 'err_floor', 1, 'read_floor', 1, 'override', NaN, 'chosen', 3);
        %params.Results.read_override_UMI_denoised=3;
        thresholds = tag_collection_denoised.compute_thresholds(params, FQ);
        tag_called_allele = tag_collection_denoised.call_alleles(thresholds.chosen, aligned);
        summary = BulkExperimentReport(tag_collection_denoised, tag_denoise_map, tag_called_allele, FQ, thresholds);
    else        
        tag_collection = CBCollection.FromFQ(FQ);
        [tag_collection_denoised, tag_denoise_map] = tag_collection.denoise(ref_CBs);
        thresholds = tag_collection_denoised.compute_thresholds(params, FQ, length(ref_CBs));
        tag_called_allele = tag_collection_denoised.call_alleles([thresholds.CB.chosen, thresholds.UMI.chosen], aligned);
        summary = SCExperimentReport(tag_collection_denoised, tag_collection, tag_denoise_map, ...
                                     tag_called_allele, FQ, thresholds, ref_CBs);
    end
    
    fprintf('Saving results...%d/%d tags edited, %d alleles\n', ...
        summary.N.eventful_tags, summary.N.called_tags, size(summary.alleles, 1));
    
    % Save just summary values needed for further analysis separately, so
    % they can be opened quickly. Full results saved later, can be big and
    % clunky to reopen for one-off analysis
    if (strcmp(cfg.type, 'Bulk'))    
        save(sprintf('%s/Summary.mat', params.Results.outdir), 'summary', 'thresholds', 'params');
    else
        save(sprintf('%s/Summary.mat', params.Results.outdir), 'summary', 'thresholds', 'params', 'ref_CBs');
    end
        
    try
        save(sprintf('%s/Analysis.mat', params.Results.outdir));
    catch
        save(sprintf('%s/Analysis.mat', params.Results.outdir), '-v7.3', '-nocompression');
    end
    gzip(sprintf('%s/Analysis.mat', params.Results.outdir));
    delete(sprintf('%s/Analysis.mat', params.Results.outdir));
    
    generate_text_output(summary, params, thresholds, params.Results.outdir);
    
    fprintf('Generating allele plot\n');
    
    warning('off', 'MATLAB:hg:AutoSoftwareOpenGL');
    plot_summary(summary, params.Results.outdir);
    
    fprintf('Generating diagnostic plot\n');
    if (strcmp(cfg.type, 'Bulk'))    
        suspect_alleles = plot_diagnostic(cfg, FQ, aligned, tag_collection_denoised, tag_denoise_map, tag_called_allele, ...
                                          summary, thresholds, params.Results.outdir);
    else
        suspect_alleles = plot_diagnostic(cfg, FQ, aligned, tag_collection_denoised, tag_denoise_map, tag_called_allele, ...
                                          summary, thresholds, ref_CBs, params.Results.outdir);
    end
    warning('on', 'MATLAB:hg:AutoSoftwareOpenGL');
                                  
    generate_warnings(summary, params, suspect_alleles, params.Results.outdir);
    
    fprintf('Pipeline completed in %g seconds\n', toc);
    
    diary off;
    copyfile(get(0,'DiaryFile'), [params.Results.outdir '/Log.txt']);
    

    
end