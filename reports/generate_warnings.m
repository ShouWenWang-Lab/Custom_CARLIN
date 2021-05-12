function generate_warnings(summary, params, suspect_alleles, output_path)

    fprintf('Generating warnings file\n');
    
    assert(isa(summary, 'ExperimentReport'));
    
    filename = [tempname '.txt'];
    fid = fopen(filename, 'wt');    
    
    sc = startsWith(params.Results.cfg_type, 'sc');
    indrops = false;
    if (sc)
        tag = 'CB';
        indrops = contains(lower(params.Results.cfg_type), 'indrop');
    else
        tag = 'UMI';
    end
    
    fprintf(fid, 'OFF-TARGET AMPLIFICATION\n');
    
    ota_pct = 100-double(summary.reads.CARLIN_match)/double(summary.reads.in_fastq)*100;
    if (ota_pct > 10)
        fprintf(fid, '\nSignificant off-target amplification detected: %.0f%% of reads are not CARLIN.\n', round(ota_pct));
    else
        fprintf(fid, ['\nInsignificant off-target amplification detected. Only %.0f%% of reads are ', ...
                      'off-target.\n'], round(ota_pct));
    end
    
    if (sc)
        ref_list_issues = false;
        fprintf(fid, '\nREFERENCE_LIST\n');
        if (any(strcmp(params.UsingDefaults, 'ref_CB_file')))
            ref_list_issues = true;
            fprintf(fid, '\nCollapsing MiSeq cell barcodes against the platform''s reference list may lead to more FPs and FNs.\n');
        else
            fprintf(fid, '\nNo issues detected with reference list.\n');
        end
    end
    
    fprintf(fid, '\nFILTERING\n');
    filtering_issues = false;
    
    low_sequencing_depth = false;    
    if (sc && summary.reads.in_fastq < 500000 || ~sc && summary.reads.in_fastq < 100000)
        fprintf(fid, '\nSequencing depth insufficent. Low number of reads detected in FASTQ.\n');
        filtering_issues = true;
        low_sequencing_depth = true;    
    end
    
    pct_CARLIN_end_found = double(summary.reads.valid_2_seq)/double(summary.reads.in_fastq)*100;
    if (ota_pct <= 10)
        if (pct_CARLIN_end_found < 75)
            errstr = sprintf('The end of the CARLIN construct is detected in only %.0f%% of reads.', round(pct_CARLIN_end_found));
            if (sc)
                errstr = [errstr ' This is likely because QC degrades toward the end of the read and is discarded by Illumina.'];                
                if (indrops)
                    errstr = [errstr ' This is a known issue with InDrops.'];            
                end
            else
                errstr = [errstr ' Consider examining the logs from PEAR to see if there were issues detected in the single-end reads.'];
            end
            fprintf(fid, '\n%s\n', errstr);
            filtering_issues = true;
        else
            pct_valid_SEQ_structure = double(summary.reads.valid_SEQ_structure)/double(summary.reads.in_fastq)*100;
            if (pct_valid_SEQ_structure < 0.75)
                fprintf(fid, ['\nOnly %.0f%% of reads have the expected structure with primers flanking CARLIN.', ...
                              'See Results.txt for a more detailed breakdown of which primer may be causing issues.\n'], ...
                              round(pct_valid_SEQ_structure));
                filtering_issues = true;
            end
        end
    end
    
    pct_usable_prov = double(summary.reads.valid_provenance_structure)/double(summary.reads.in_fastq)*100;
    low_usable_prov = pct_usable_prov < 75;
    if (low_usable_prov)
        fprintf(fid, ['\nOnly %.0f%% of reads have usable provenance information (CB or UMI). ' ...
                      'See Results.txt for a more detailed breakdown of QC issues.\n'], round(pct_usable_prov));
        filtering_issues = true;
    end
    
    if (~filtering_issues)
        fprintf(fid, '\nNo issues detected at filtering step.\n');
    end
    
    fprintf(fid, '\nANALYSIS\n');
    analysis_issues = false;
  
    reads_per_edited_tag = double(summary.reads.eventful_tags_total)/double(summary.N.eventful_tags);
    reads_per_template_tag = double(summary.reads.called_tags_total-summary.reads.eventful_tags_total)/...
                             double(summary.N.called_tags-summary.N.eventful_tags);
    pref_amp_detected = reads_per_edited_tag / reads_per_template_tag > 2;
    
    if (pref_amp_detected)
        fprintf(fid, '\nPreferential amplification of edited alleles detected. (Reads per edited %s / reads per unedited %s) > 2.\n', tag, tag);
        analysis_issues = true;
    end
        
    if (summary.N.common_tags < 500)
        prestr = sprintf('Number of common %ss is low (%d). ', tag, summary.N.common_tags);
    else
        prestr = '';
    end
    pct_common_reads = double(summary.reads.common_tags)/double(summary.reads.valid_lines)*100;
    if (pct_common_reads < 70)
        prestr = [prestr sprintf('Common %ss account for only %.0f%% of reads. ', tag, round(pct_common_reads))];
    end
    
    if (~isempty(prestr))
        analysis_issues = true;
        reason_known = false;
        poststr = ['This is likely due to:'];
        if (sc && ref_list_issues)
            reason_known = true;
            poststr = sprintf('%s\n - reads dedicated to FP CBs owing to the inflated reference list used', poststr);
        end
        if (low_sequencing_depth)
            reason_known = true;
            poststr = sprintf('%s\n - low sequencing depth so that many %ss do not exceed the required threshold', poststr, tag);
        end
        if (pref_amp_detected)
            reason_known = true;
            poststr = sprintf('%s\n - preferential amplification of edited transcripts causing many reads to be allocated to a few %ss', poststr, tag);
        end
        if (low_usable_prov)
            reason_known = true;
            poststr = sprintf(['%s\n - QC issues with CB/UMI which cause reads to be dedicated to ' ...
                               'spurious %ss that persist after filtering'], poststr, tag);
        end
        if (~reason_known)
            poststr = 'Reason unknown.';
        end
        fprintf(fid, '\n%s%s\n', prestr, poststr);
    end        
    
    pct_callable_tags = double(summary.N.called_tags)/double(summary.N.common_tags)*100;
    pct_callable_reads = double(summary.reads.called_tags_allele)/double(summary.reads.called_tags_total) * 100;    
    low_call_yield = pct_callable_tags < 70 || pct_callable_reads < 70;
    errstr = '';
    if (pct_callable_tags < 70)
        errstr = sprintf('Only %.0f%% of common %ss could successfully call an allele. ', round(pct_callable_tags), tag);
    end
    if (pct_callable_reads < 70)
        errstr = [errstr sprintf('Only %0.f%% of reads in callable %ss folded into consensus for allele. ', round(pct_callable_reads), tag)];
    end
    
    if (low_call_yield)
        analysis_issues = true;
        fprintf(fid, ['\n%sThis suggests an issue with library preparation, such as low primer specificity ', ...
                      'or cross-talk during amplification.\n'], errstr);
    end
    
    if (~analysis_issues)
        fprintf(fid, '\nNo issues detected during analysis.\n');
    end
        
    fprintf(fid, '\nRESULTS\n');
    results_issues = false;
      
    pct_edited = double(summary.N.eventful_tags)/double(summary.N.called_tags)*100;
    CP = mean(cellfun(@(x) CARLIN_def.getInstance.N.segments-length(Mutation.find_modified_sites(x)), summary.alleles));
    
    if (pct_edited < 30)
        fprintf(fid, '\nLow +Dox induction detected. Only %.0f%% of %ss reported an edited CARLIN allele.\n', round(pct_edited), tag);
        results_issues = true;
    else    
        template_mask = strcmp(cellfun(@(x) degap(x.get_seq), summary.alleles, 'un', false), ...
                               CARLIN_def.getInstance.seq.CARLIN);
        prestr = [];
        if (~any(template_mask))
            prestr = sprintf('Template allele never reported by any %s.', tag);
        elseif (find(template_mask)~=1)
            prestr = 'Template allele was not the most common allele.';
        end
        if (~isempty(prestr))
            results_issues = true;
            if (low_call_yield && pref_amp_detected)
                fprintf(fid, ['\n%s This is likely because preferential amplification ' ...
                              'of short sequences and difficulty in establishing consensus for many %ss causes '...
                              'the template to be under-represented in reads\n'], prestr, tag);
            elseif (low_call_yield)
                fprintf(fid, ['\n%s This is likely because difficulty in establishing allele-calling consensus for many %ss '...
                              'penalizes longer sequences\n'], prestr, tag);
            elseif (pref_amp_detected)
                fprintf(fid, ['\n%s This is likely because preferential amplification ' ...
                              'of short sequences caused the template to be under-represented after PCR\n'], prestr);
            elseif (CP < 3)
                fprintf(fid, ['\n%s Since significant preferential amplification ' ...
                              'of short sequences was not detected and most common %ss called alleles successfully, ' ...
                              'this may indicate over-exposure to +Dox as most alleles are highly degraded\n'], prestr, tag);
            else
                fprintf(fid, ['\n%s This is not because of preferential amplification, low yield in calling alleles among common %ss ' ...
                              'or over-exposure to +Dox.\n'], prestr, tag);
            end
        end
    end
    
    if (~isempty(suspect_alleles))
        results_issues = true;        
        if (sc)
            errstr = sprintf(['\nFor allele(s) %s, >10%% of CB halves or UMIs differ pairwise by only 1bp. ' ...
                     'This may indicate that that the cell count for these alleles is artificially inflated.'], num2str(suspect_alleles));
        else
            errstr = sprintf(['\nFor allele(s) %s, >10%% of UMIs differ pairwise by only 1bp. ' ...
                     'This may indicate that that the UMI count for these alleles is artificially inflated.'], num2str(suspect_alleles));
        end
        if (low_usable_prov && length(suspect_alleles)/length(summary.alleles)*100>5)            
            errstr = [errstr ' This is likely because of QC issues on the CB/UMI reads.'];
        end
        if (sc && indrops)
            errstr = [errstr ' This is a known issue with InDrops.'];
        end
        fprintf(fid, '%s\n', errstr);
    end
    
    if (~results_issues)
        fprintf(fid, '\nNo issues detected in results.\n');
    end
    
    fclose(fid);
    copyfile(filename, [output_path '/Warnings.txt']);
    
end