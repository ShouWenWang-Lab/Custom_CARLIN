function tests = CARLIN_pipeline_test
    tests = functiontests(localfunctions);
end

function test_CARLIN_pipeline_DNA(testCase)
    check_CARLIN_pipeline_bulk(testCase, 'DNA')
end

function test_CARLIN_pipeline_RNA(testCase)
    check_CARLIN_pipeline_bulk(testCase, 'RNA')
end

function test_CARLIN_pipeline_InDropsV2(testCase)    
    check_CARLIN_pipeline_InDrops(testCase, 2)
end

function test_CARLIN_pipeline_InDropsV3(testCase)    
    check_CARLIN_pipeline_InDrops(testCase, 3)
end

function test_CARLIN_pipeline_10xV2(testCase)
    check_CARLIN_pipeline_10x(testCase, 2)
end

function test_CARLIN_pipeline_10xV3(testCase)
    check_CARLIN_pipeline_10x(testCase, 3)
end

function check_CARLIN_pipeline_bulk(testCase, bulk_type)    
    [folder, ~, ~] = fileparts(mfilename('fullpath'));  
    cfg_type = sprintf('Bulk%s', bulk_type);
    out_path = sprintf('%s/Output/%s', folder, cfg_type);
    if (exist(out_path, 'dir')==7)
        rmdir(out_path, 's');
    end
    fastq_file = sprintf('%s/data/%s.fastq.gz', folder, cfg_type);    
    verifyWarningFree(testCase, @() analyze_CARLIN(fastq_file, cfg_type, out_path));
    check_all_files(testCase, out_path);
    gunzip(sprintf('%s/Analysis.mat.gz', out_path));
    load(sprintf('%s/Analysis.mat', out_path));
    check_bulk_FQ(testCase, FQ, fastq_file, cfg);
    check_bulk_collection(testCase, FQ, tag_collection, containers.Map(keys(tag_denoise_map), keys(tag_denoise_map)));
    check_bulk_collection(testCase, FQ, tag_collection_denoised, tag_denoise_map);
    check_plaintext_results(testCase, summary, out_path);
end

function check_CARLIN_pipeline_InDrops(testCase, version)
    [folder, ~, ~] = fileparts(mfilename('fullpath'));
    cfg_type = sprintf('scInDropsV%d', version);    
    out_path = sprintf('%s/Output/%s', folder, cfg_type);
    if (exist(out_path, 'dir')==7)
        rmdir(out_path, 's');
    end
    fastq_file = sprintf('%s/data/InDropsV%d.fastq.gz', folder, version);
    verifyWarningFree(testCase, @() analyze_CARLIN(fastq_file, cfg_type, out_path));
    check_all_files(testCase, out_path);
    gunzip(sprintf('%s/Analysis.mat.gz', out_path));
    load(sprintf('%s/Analysis.mat', out_path), 'FQ', 'cfg');
    check_sc_FQ(testCase, FQ, fastq_file, cfg);
end

function check_CARLIN_pipeline_10x(testCase, version)
    [folder, ~, ~] = fileparts(mfilename('fullpath'));
    cfg_type = sprintf('sc10xV%d', version);
    out_path = sprintf('%s/Output/%s', folder, cfg_type);
    if (exist(out_path, 'dir')==7)
        rmdir(out_path, 's');
    end
    fastq_file = {sprintf('%s/data/10xV%d_R1.fastq.gz', folder, version), ...
                  sprintf('%s/data/10xV%d_R2.fastq.gz', folder, version)};
    verifyWarningFree(testCase, @() analyze_CARLIN(fastq_file, cfg_type, out_path));
    gunzip(sprintf('%s/Analysis.mat.gz', out_path));
    check_all_files(testCase, out_path);
    load(sprintf('%s/Analysis.mat', out_path), 'FQ', 'cfg');
    check_sc_FQ(testCase, FQ, fastq_file, cfg);
end

function check_bulk_FQ(testCase, FQ, fastq_file, cfg)
    [SEQ, read_SEQ, QC] = BulkFastQData.parse_bulk_fastq(fastq_file, cfg);    
    [UMI_check, QC_check] = FastQData.orient_reads(cfg, FQ.UMI, FQ.QC);
    verifyTrue(testCase, all(cellfun(@(x,y) contains(x, char(y)), QC, QC_check)));
    verifyTrue(testCase, all(cellfun(@(x,y) contains(x, y), SEQ(read_SEQ), UMI_check(FQ.read_UMI))));
    verifyTrue(testCase, all(cellfun(@(x,y) contains(x, int2nt(y)), SEQ(read_SEQ), FQ.SEQ_raw(FQ.read_SEQ_raw))));
    verifyTrue(testCase, all(cellfun(@(x,y,s,e1,e2) isequal(x(s:min(e1,e2)), y), ...
                            FQ.SEQ_raw(FQ.read_SEQ_raw(FQ.masks.valid_lines)), ...
                            FQ.SEQ_trimmed(FQ.read_SEQ_trimmed(FQ.masks.valid_lines)), ...
                            num2cell(FQ.trim_loc.head_after_trim_5_primer(FQ.read_SEQ_raw(FQ.masks.valid_lines))), ...
                            num2cell(FQ.trim_loc.tail_after_trim_3_primer(FQ.read_SEQ_raw(FQ.masks.valid_lines))), ...
                            num2cell(FQ.trim_loc.tail_after_trim_2_seq(FQ.read_SEQ_raw(FQ.masks.valid_lines))))));
    verifyEqual(testCase, FQ.SEQ_trimmed(FQ.read_SEQ_trimmed(FQ.masks.valid_lines)), ...
                          FQ.SEQ_valid  (FQ.read_SEQ_valid  (FQ.masks.valid_lines)));    
end

function check_sc_FQ(testCase, FQ, fastq_file, cfg)
    if (strcmp(cfg.SC.Platform, 'InDrops'))
        H   = cell(FQ.Nreads,1);
        SEQ = cell(FQ.Nreads,1);
        [fastq_file, ext] = FastQData.maybe_unzip(fastq_file);
        [H(:,1), SEQ(:,1), ~] = fastqread(fastq_file);
        FastQData.maybe_clear_unzipped(fastq_file, ext);
        H = H(FQ.masks.valid_lines);
        SEQ = SEQ(FQ.masks.valid_lines);
        [CB, UMI, QC] = SCFastQData.parse_indrops_provenance(H); 
        clear H;        
    elseif(strcmp(cfg.SC.Platform, '10x'))
        CB = cell(FQ.Nreads,1);
        QC = cell(FQ.Nreads,1);
        SEQ = cell(FQ.Nreads,1);
        [fastq_file, ext] = FastQData.maybe_unzip(fastq_file);
        [~, CB(:,1), QC(:,1)] = fastqread(fastq_file{1});
        [~, SEQ(:,1), ~] = fastqread(fastq_file{2});
        FastQData.maybe_clear_unzipped(fastq_file, ext);
        CB = CB(FQ.masks.valid_lines);
        QC = QC(FQ.masks.valid_lines);
        SEQ = SEQ(FQ.masks.valid_lines);
        [CB, UMI, QC] = SCFastQData.parse_10x_provenance(CB, QC, cfg);
    end
    verifyEqual(testCase, FQ.CB(FQ.read_CB(FQ.masks.valid_lines)), CB);
    verifyEqual(testCase, FQ.UMI(FQ.read_UMI(FQ.masks.valid_lines)), UMI);
    verifyEqual(testCase, cellfun(@(x) int2nt(x), FQ.SEQ_raw(FQ.read_SEQ_raw(FQ.masks.valid_lines)), 'un', false), SEQ);
    verifyEqual(testCase, FQ.QC(FQ.masks.valid_lines), cellfun(@uint8, QC, 'un', false));
    verifyTrue(testCase, all(cellfun(@(x,y) contains(int2nt(x), int2nt(y)), ...
                            FQ.SEQ_raw(FQ.read_SEQ_raw(FQ.masks.valid_lines)), ...
                            FQ.SEQ_trimmed(FQ.read_SEQ_trimmed(FQ.masks.valid_lines)))));
    verifyEqual(testCase, FQ.SEQ_trimmed(FQ.read_SEQ_trimmed(FQ.masks.valid_lines)), ...
                          FQ.SEQ_valid  (FQ.read_SEQ_valid  (FQ.masks.valid_lines)));    

end

function check_bulk_collection(testCase, FQ, tag_collection, UMI_map)
    k = keys(UMI_map);
    [~, which_UMI] = ismember(k, FQ.UMI);
    v = values(UMI_map);
    [v, ~, ind] = unique(v);
    [~, where] = ismember({tag_collection.UMIs.UMI}', v);
    source = arrayfun(@(i) which_UMI(ind==i), where, 'un', false);    
    
    which_SEQ = {tag_collection.UMIs.SEQ_ind}';
    source = repelem(source, cellfun(@length, which_SEQ));
    which_SEQ = num2cell(vertcat(which_SEQ{:}));
    SEQ_weight = num2cell(vertcat(tag_collection.UMIs.SEQ_weight));
    
    verifyTrue(testCase, all(cellfun(@(u,s,w) ...
               sum(ismember(FQ.read_UMI(FQ.masks.valid_lines), u) & FQ.read_SEQ_valid(FQ.masks.valid_lines)==s) == w, ...
               source, which_SEQ, SEQ_weight)));
end

function check_plaintext_results(testCase, summary, out_path)
    mut_list = cellfun(@Mutation.identify_Cas9_events, summary.alleles, 'un', false);
    alleles_rec = cellfun(@Mutation.apply_mutations, mut_list, 'un', false);
    alleles = cellfun(@Mutation.apply_mutations, Mutation.FromFile([out_path '/AlleleAnnotations.txt']), 'un', false);
    verifyEqual(testCase, alleles, alleles_rec);
    colonies = cellfun(@(x) strsplit(x, ',')', splitlines(fileread([out_path '/AlleleColonies.txt'])), 'un', false);
    verifyEqual(testCase, colonies, summary.allele_colony)
end

function check_all_files(testCase, outpath)
    verifyEqual(testCase, exist(sprintf('%s/Summary.mat',           outpath), 'file'), 2);
    verifyEqual(testCase, exist(sprintf('%s/Analysis.mat.gz',       outpath), 'file'), 2);
    verifyEqual(testCase, exist(sprintf('%s/Alleles.png',           outpath), 'file'), 2);
    verifyEqual(testCase, exist(sprintf('%s/Diagnostic.png',        outpath), 'file'), 2);
    verifyEqual(testCase, exist(sprintf('%s/AlleleAnnotations.txt', outpath), 'file'), 2);
    verifyEqual(testCase, exist(sprintf('%s/AlleleColonies.txt',    outpath), 'file'), 2);
    verifyEqual(testCase, exist(sprintf('%s/Log.txt',               outpath), 'file'), 2);
    verifyEqual(testCase, exist(sprintf('%s/Results.txt',           outpath), 'file'), 2);
    verifyEqual(testCase, exist(sprintf('%s/Warnings.txt',          outpath), 'file'), 2);
end

