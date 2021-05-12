function suspect_alleles = plot_diagnostic(cfg, FQ, aligned, tag_collection_denoised, tag_denoise_map, ...
                                           tag_called_allele, summary, thresholds, varargin)
    
    [~, ~, event_ind] = unique(cellfun(@(x) x.get_event_structure, aligned.aligned_SEQ, 'un', false));
    N = length(FQ.masks.valid_lines);
    
    dat = struct('sequence' , FQ.read_SEQ_valid(FQ.masks.valid_lines), ...
                 'sanitized', uint32(aligned.alignment_map(FQ.read_SEQ_valid(FQ.masks.valid_lines))), ...
                 'event'    , uint32(event_ind(aligned.alignment_map(FQ.read_SEQ_valid(FQ.masks.valid_lines)))), ...
                 'allele'   , uint32(zeros(N,1)), ...
                 'allele_eventual', uint32(zeros(N,1)));
    
    if (strcmp(cfg.type, 'Bulk'))
        
        [is, k] = ismember(keys  (tag_denoise_map), FQ.UMI);
        assert(all(is));
        [is, v] = ismember(values(tag_denoise_map), FQ.UMI);
        assert(all(is));
        
        dat.CB_raw = FQ.read_UMI(FQ.masks.valid_lines);
        dat.CB = dat.CB_raw;
        [is, where] = ismember(dat.CB_raw, k);
        dat.CB(is) = v(where(is));
        
        if (size(summary.alleles,1)>0)
        
            [is, which_UMI] = ismember(vertcat(summary.allele_colony{:}), FQ.UMI);
            assert(all(is));
            which_UMI = mat2cell(which_UMI, cellfun(@length, summary.allele_colony));

            [is, where_in_collection] = ismember(vertcat(summary.allele_colony{:}), {tag_collection_denoised.UMIs.UMI}');
            assert(all(is));
            where_in_collection = mat2cell(where_in_collection, cellfun(@length, summary.allele_colony));

            [is, where_in_called] = ismember(vertcat(summary.allele_colony{:}), {tag_called_allele.UMI}');
            assert(all(is));
            where_in_called = mat2cell(where_in_called, cellfun(@length, summary.allele_colony));

            for i = 1:size(summary.alleles,1)
                for j = 1:size(summary.allele_colony{i})                
                    which_SEQ = tag_collection_denoised.UMIs(where_in_collection{i}(j)).SEQ_ind(...
                                tag_called_allele(where_in_called{i}(j)).constituents);
                    how_much = tag_collection_denoised.UMIs(where_in_collection{i}(j)).SEQ_weight(...
                               tag_called_allele(where_in_called{i}(j)).constituents);
                    dat.allele((dat.CB==which_UMI{i}(j)) & ismember(dat.sequence,which_SEQ)) = i;
                    dat.allele_eventual(dat.CB==which_UMI{i}(j)) = i;
                    assert(sum((dat.CB==which_UMI{i}(j)) & ismember(dat.sequence,which_SEQ))==sum(how_much));
                end
            end
        end
        assert(all(dat.allele_eventual(dat.allele>0)==dat.allele(dat.allele>0)));
        
    else    
        
        assert(~isempty(varargin));
        
        ref_list = varargin{1};
        
        [is, k] = ismember(keys  (tag_denoise_map.CB), FQ.CB);
        assert(all(is));
        [is, v] = ismember(values(tag_denoise_map.CB), ref_list);
        assert(all(is));
        dat.CB_raw = FQ.read_CB(FQ.masks.valid_lines);
        dat.CB     = zeros(size(dat.CB_raw));
        [is, where] = ismember(dat.CB_raw, k);
        dat.CB(is) = v(where(is));
        v = unique(v);
        
        dat.UMI_raw = FQ.read_UMI(FQ.masks.valid_lines);
        dat.UMI     = zeros(size(dat.UMI_raw));
        
        kUMI = cellfun(@(x) keys(tag_denoise_map.UMI(x))', ref_list(v), 'un', false);
        vUMI = cellfun(@(x) values(tag_denoise_map.UMI(x))', ref_list(v), 'un', false);
        nUMI = cellfun(@length, vUMI);
        kUMI = vertcat(kUMI{:});
        vUMI = vertcat(vUMI{:});
        [is, kUMI] = ismember(kUMI, FQ.UMI);
        assert(all(is));
        [is, vUMI] = ismember(vUMI, FQ.UMI);
        assert(all(is));
        kUMI = mat2cell(kUMI, nUMI);
        vUMI = mat2cell(vUMI, nUMI);
        for i = 1:length(v)
            [is, where] = ismember(dat.UMI_raw(dat.CB==v(i)), kUMI{i});
            assert(all(is));
            dat.UMI(dat.CB==v(i)) = vUMI{i}(where(is));
        end
        assert(sum(dat.UMI==0)==sum(dat.CB==0));
        
        if (size(summary.alleles,1)>0)
        
            [is, which_CB] = ismember(vertcat(summary.allele_colony{:}), ref_list);
            assert(all(is));
            which_CB = mat2cell(which_CB, cellfun(@length, summary.allele_colony));

            [is, where_in_collection] = ismember(vertcat(summary.allele_colony{:}), {tag_collection_denoised.CBs.CB}');
            assert(all(is));
            where_in_collection = mat2cell(where_in_collection, cellfun(@length, summary.allele_colony));

            [is, where_in_called] = ismember(vertcat(summary.allele_colony{:}), {tag_called_allele.CB}');
            assert(all(is));
            where_in_called = mat2cell(where_in_called, cellfun(@length, summary.allele_colony));

            for i = 1:size(summary.alleles,1)
                for j = 1:size(summary.allele_colony{i})
                    supporting_UMIs = tag_called_allele(where_in_called{i}(j)).constituents;
                    which_UMIs = tag_collection_denoised.CBs(where_in_collection{i}(j)).UMIs(supporting_UMIs);
                    CB_mask = (dat.CB==which_CB{i}(j));
                    for k = 1:size(which_UMIs,1)
                        which_SEQ = tag_called_allele(where_in_called{i}(j)).umi_call_result(supporting_UMIs(k)).constituents;
                        how_much  = which_UMIs(k).SEQ_weight(which_SEQ);
                        which_SEQ = which_UMIs(k).SEQ_ind(which_SEQ);
                        set_mask = CB_mask & ismember(dat.sequence,which_SEQ) & (dat.UMI==find(strcmp(FQ.UMI, which_UMIs(k).UMI)));
                        assert(sum(set_mask)==sum(how_much));
                        dat.allele(set_mask) = i;                    
                    end
                    dat.allele_eventual(CB_mask) = i;
                end
            end
        end
        assert(all(dat.allele_eventual(dat.allele>0)==dat.allele(dat.allele>0)));
        discard_mask = dat.CB==0;
        assert(sum(~discard_mask) == summary.reads.matched_tags);
        dat = structfun(@(x) x(~discard_mask), dat, 'un', false);        
    end
    
    [r, ~, v] = find(accumarray(nonzeros(dat.CB), 1));    
    [v, idx] = sort(v, 'descend');
    [~, dat.CB] = ismember(dat.CB, r(idx));
    
    if (strcmp(cfg.type, 'Bulk'))
        cutoff = find(v>=thresholds.chosen, 1, 'last');
    else
        cutoff = find(v>=thresholds.CB.chosen, 1, 'last');
    end
    if (isempty(cutoff))
        cutoff = 1;
    else
        assert(cutoff == summary.N.common_tags);
    end
   
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    set(fig,'defaultAxesColorOrder',[[0, 0, 0]; [0, 0, 0]]);
    
    Nr = 4;
    Nc = 5;
    
    if (strcmp(cfg.type, 'Bulk'))        
        sp(1,1) = reads_by_filter    (Nr, Nc, 1, summary, 'UMI');        
        sp(1,2) = plot_read_distribution(Nr, Nc, 2, tag_collection_denoised, thresholds, summary, 'UMI');
        sp(1,3) = length_distribution(Nr, Nc, 3, FQ, summary, 'UMI');
        sp(1,4) = UMI_composition_all(Nr, Nc, 4, FQ, summary);                           
        sp(1,5) = UMI_composition_by_allele(Nr, Nc, 5, summary);
        sp(2,1) = reads_by_tag       (Nr, Nc, [ 6: 9], dat, cutoff, summary, 'UMI');
        sp(3,1) = distance_by_tag    (Nr, Nc, [11:14], dat, cutoff, FQ, aligned, summary);
        sp(4,1) = variety_by_tag     (Nr, Nc, [16:19], dat, cutoff, 'UMI');        
        sp(2,2) = reads_by_allele    (Nr, Nc, 10, dat, summary);
        [sp(3,2), suspect_alleles] = distance_by_allele (Nr, Nc, 15, dat, FQ, aligned, summary, 'UMI');
        sp(4,2) = variety_by_allele  (Nr, Nc, 20, dat, 'UMI');                     
    else
        sp(1,1) = reads_by_filter(Nr, Nc, 1, summary, 'CB');
        sp(1,2) = plot_read_distribution(Nr, Nc, 2, tag_collection_denoised, thresholds, summary, 'CB');
        sp(1,3) = length_distribution(Nr, Nc, 3, FQ, summary, 'CB');         
        sp(1,4) = CB_UMI_composition_all(Nr, Nc, 4, FQ, summary, dat);
        sp(1,5) = CB_UMI_composition_by_allele(Nr, Nc, 5, summary, FQ, dat);
        sp(2,1) = reads_by_tag       (Nr, Nc, [ 6: 9], dat, cutoff, summary, 'CB');
        sp(3,1) = distance_by_tag    (Nr, Nc, [11:14], dat, cutoff, FQ, aligned, summary);
        sp(4,1) = variety_by_tag     (Nr, Nc, [16:19], dat, cutoff, 'CB');
        sp(2,2) = reads_by_allele    (Nr, Nc, 10, dat, summary);
        [sp(3,2), suspect_alleles] = distance_by_allele (Nr, Nc, 15, dat, FQ, aligned, summary, 'CB');
        sp(4,2) = variety_by_allele  (Nr, Nc, 20, dat, 'CB');        
    end
    
    warning('off', 'MATLAB:linkaxes:RequireDataAxes');
    linkaxes(sp(2:4,1), 'x');
    linkaxes(sp(2:4,2), 'x');
    linkaxes(sp(3,:), 'y');
    
    position_subplots(cutoff, size(summary.alleles,1), sp);
    
    if ((strcmp(cfg.type, 'Bulk') && length(varargin)==1) || ...
        (strcmp(cfg.type, 'SC') && length(varargin)==2))
        outdir = varargin{end};    
        if (~exist(outdir, 'dir'))
            mkdir(outdir);
        end
        print(fig, sprintf('%s/Diagnostic.png', outdir), '-dpng','-r0');
        close;
    end
    
end

function sp = plot_read_distribution(Nr, Nc, which_sp, CB_collection_denoised, thresholds, summary, type)

    sp = subplot(Nr, Nc, which_sp);
    if (strcmp(type, 'UMI'))
        freq = CB_collection_denoised.get_UMI_freq();        
        plot(log10(1:length(freq)), log10(freq), 'Color', 'black'); hold on;        
        x_cut = log10(summary.N.common_tags);
        y_cut = log10(thresholds.chosen);
        line([x_cut, x_cut], [0 y_cut], 'LineStyle', ':');
        line([0 x_cut], [y_cut y_cut], 'LineStyle', ':');        
        hold off;
        lgd = legend({'UMI'}, 'Location', 'Southwest');
        lgd.FontSize = 6;
    else
        CB_freq = CB_collection_denoised.get_CB_freq();
        UMI_freq = CB_collection_denoised.get_UMI_freq();        
        plot(log10(1:length(CB_freq)), log10(CB_freq), 'Color', 'black'); hold on;
        plot(log10(1:length(UMI_freq)), log10(UMI_freq), 'Color', 'blue');        
        x_cut = log10(summary.N.common_tags);
        y_cut = log10(thresholds.CB.chosen);
        line([x_cut, x_cut], [0 y_cut], 'Color', 'black', 'LineStyle', ':');
        line([0 x_cut], [y_cut y_cut], 'Color', 'black', 'LineStyle', ':');
        x_cut = log10(find(UMI_freq>=thresholds.UMI.chosen, 1, 'last'));
        y_cut = log10(thresholds.UMI.chosen);
        line([x_cut, x_cut], [0 y_cut], 'Color', 'blue', 'LineStyle', ':');
        line([0 x_cut], [y_cut y_cut], 'Color', 'blue', 'LineStyle', ':');
        hold off;
        lgd = legend({'CB'; 'UMI'}, 'Location', 'Southwest');
        lgd.FontSize = 6;
    end
        
    xlabel('log10 Rank');
    ylabel('log10 Reads');
    set(gca,'XAxisLocation', 'top');
    yrule = get(gca, 'YAxis');
    yrule.FontSize = 6;        
    axis tight; box on;
    
end

function sp = reads_by_filter(Nr, Nc, which_sp, summary, tag)

    fn = fieldnames(summary.reads);
    num_reads = cellfun(@(x) getfield(summary.reads, x), fn, 'un', false);
    num_reads = vertcat(num_reads{:});
        
    sp = subplot(Nr, Nc, which_sp);
    
    barh(num_reads, 'blue');
    xlabel('Reads');
    yticks(1:length(num_reads));
    yticklabels(strrep(fn, 'tag', tag));
    set(gca,'TickLabelInterpreter','none', 'XAxisLocation', 'top', 'YDir', 'reverse');
    yrule = get(gca, 'YAxis');
    yrule.FontSize = 6;    
    axis tight; box on;
end

function sp = length_distribution(Nr, Nc, which_sp, FQ, summary, type)

    L_raw_read     = accumarray(cellfun(@length, FQ.SEQ_raw)+1, accumarray(nonzeros(FQ.read_SEQ_raw), 1));
    if (strcmp(type, 'UMI'))
        raw_read_deduped = unique([FQ.read_SEQ_raw FQ.read_UMI], 'rows');        
    else
        raw_read_deduped = unique([FQ.read_SEQ_raw FQ.read_CB FQ.read_UMI], 'rows');        
    end
    L_raw_norm     = accumarray(cellfun(@length, FQ.SEQ_raw)+1, accumarray(nonzeros(raw_read_deduped(:,1)), 1));

    L_trimmed_read = accumarray(cellfun(@length, FQ.SEQ_trimmed)+1, accumarray(nonzeros(FQ.read_SEQ_trimmed), 1));
    if (strcmp(type, 'UMI'))
        trimmed_read_deduped = unique([FQ.read_SEQ_trimmed FQ.read_UMI], 'rows');        
    else
        trimmed_read_deduped = unique([FQ.read_SEQ_trimmed FQ.read_CB FQ.read_UMI], 'rows');        
    end
    L_trimmed_norm = accumarray(cellfun(@length, FQ.SEQ_trimmed)+1, accumarray(nonzeros(trimmed_read_deduped(:,1)), 1));
    
    L_valid_read   = accumarray(cellfun(@length, FQ.SEQ_valid)+1, accumarray(nonzeros(FQ.read_SEQ_valid), 1));
    if (strcmp(type, 'UMI'))
        valid_read_deduped = unique([FQ.read_SEQ_valid FQ.read_UMI], 'rows');        
    else
        valid_read_deduped = unique([FQ.read_SEQ_valid FQ.read_CB FQ.read_UMI], 'rows');        
    end
    L_valid_norm = accumarray(cellfun(@length, FQ.SEQ_valid)+1, accumarray(nonzeros(valid_read_deduped(:,1)), 1));
    
    if (~isempty(summary.alleles))
        L_cell_norm  = accumarray(cellfun(@(x) length(degap(x.get_seq)), summary.alleles)+1, summary.allele_freqs);
    else
        L_cell_norm = [];
    end
    
    L_raw_read     = L_raw_read     / sum(L_raw_read);
    L_raw_norm     = L_raw_norm     / sum(L_raw_norm);
    L_trimmed_read = L_trimmed_read / sum(L_trimmed_read);
    L_trimmed_norm = L_trimmed_norm / sum(L_trimmed_norm);
    L_valid_read   = L_valid_read   / sum(L_valid_read);
    L_valid_norm   = L_valid_norm   / sum(L_valid_norm);
    L_cell_norm    = L_cell_norm    / sum(L_cell_norm);
    
    max_len = max([length(L_raw_read), ...
                   length(L_raw_norm), ...
                   length(L_trimmed_read), ... 
                   length(L_trimmed_norm), ... 
                   length(L_valid_read), ...
                   length(L_valid_norm), ...
                   length(L_cell_norm)]);
               
    L_raw_read     = [L_raw_read;     zeros(max_len-length(L_raw_read),    1)];
    L_raw_norm     = [L_raw_norm;     zeros(max_len-length(L_raw_norm),    1)];
    L_trimmed_read = [L_trimmed_read; zeros(max_len-length(L_trimmed_read),1)];
    L_trimmed_norm = [L_trimmed_norm; zeros(max_len-length(L_trimmed_norm),1)];
    L_valid_read   = [L_valid_read;   zeros(max_len-length(L_valid_read),  1)];
    L_valid_norm   = [L_valid_norm;   zeros(max_len-length(L_valid_norm),  1)];
    L_cell_norm    = [L_cell_norm;    zeros(max_len-length(L_cell_norm),   1)];
    
    sp = subplot(Nr, Nc, which_sp);
    
    patch([0 0 max_len-1 max_len-1], [3, 4, 4, 3], 'black', 'FaceAlpha', 0.2); hold on;
    patch([0 0 max_len-1 max_len-1], [1, 2, 2, 1], 'black', 'FaceAlpha', 0.2);
    
    plot(0:max_len-1, L_raw_read+3,     'Color', 'black');
    plot(0:max_len-1, L_raw_norm+3,     'Color', 'black', 'LineStyle', ':');
    plot(0:max_len-1, L_trimmed_read+2, 'Color', 'blue');
    plot(0:max_len-1, L_trimmed_norm+2, 'Color', 'blue', 'LineStyle', ':');
    plot(0:max_len-1, L_valid_read+1,   'Color', 'green');
    plot(0:max_len-1, L_valid_norm+1,   'Color', 'green', 'LineStyle', ':');
    plot(0:max_len-1, L_cell_norm,      'Color', 'red');  hold off;
    
    xlabel('Length');    
    xlim([0, max_len-1]);
    
    set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Ytick', [0.5:1:3.5], ...
        'Yticklabels', {'Cell'; 'Valid'; 'Trimmed'; 'FASTQ'}); 
    yrule = get(gca, 'YAxis');
    yrule.FontSize = 6;
    ytickangle(90);
    axis tight; box on;
end

function sp = UMI_composition_all(Nr, Nc, which_sp, FQ, summary)

    selector = FQ.masks.UMI_correct_length;
    hm_all   = seqprofile(FQ.UMI(FQ.read_UMI(selector)), 'Alphabet', 'NT');
    selector = intersect(selector, FQ.masks.valid_SEQ_structure);
    hm_trimmed   = seqprofile(FQ.UMI(FQ.read_UMI(selector)), 'Alphabet', 'NT');
    hm_valid = seqprofile(FQ.UMI(FQ.read_UMI(FQ.masks.valid_lines)), 'Alphabet', 'NT');
    if (~isempty(summary.allele_colony))
        hm_cell  = seqprofile(vertcat(summary.allele_colony{:}), 'Alphabet', 'NT');
    else
        hm_cell  = zeros(4, size(hm_all,2));
    end

    hm       = cumsum([hm_all; hm_trimmed; hm_valid; hm_cell], 1, 'reverse');

    sp = subplot(Nr,Nc,which_sp);
    hold on;
    for i = 1:4
        bar(hm((i-1)*4+1,:), 'blue',   'BarWidth', 0.99);
        bar(hm((i-1)*4+2,:), 'green',  'BarWidth', 0.99);
        bar(hm((i-1)*4+3,:), 'yellow', 'BarWidth', 0.99);
        bar(hm((i-1)*4+4,:), 'red',    'BarWidth', 0.99);
    end
    hold off;
    
    set(gca, 'Xtick', [], 'Ytick', []);
    axis tight; box on;
    
end

function sp = CB_UMI_composition_all(Nr, Nc, which_sp, FQ, summary, dat)
   

%     [CB1, CB2] = cellfun(@(x) deal(x(1:end-8), x(end-7:end)), FQ.CB(FQ.read_CB), 'un', false);
    [CB1, CB2] = cellfun(@(x) deal(x(1:end-8), x(end-7:end)), FQ.CB, 'un', false);
    CB1 = strrep(CB1, 'N', '-');
    CB2 = strrep(CB2, 'N', '-');    
    maxlen = max(cellfun(@length, CB1));
    CB1 = pad(CB1, maxlen, 'left', '-');
    
    CB1 = CB1(FQ.read_CB);
    CB2 = CB2(FQ.read_CB);
    
    UMI = strrep(FQ.UMI, 'N', '-');
    UMI = UMI(FQ.read_UMI);

    selector = intersect(FQ.masks.CB_correct_length, FQ.masks.UMI_correct_length);
    hm_all   = [seqprofile(CB1(selector), 'Alphabet', 'NT', 'Gaps', 'all'), ...
                seqprofile(CB2(selector), 'Alphabet', 'NT', 'Gaps', 'all'), ...
                seqprofile(UMI(selector), 'Alphabet', 'NT', 'Gaps', 'all')];
    selector = intersect(selector, FQ.masks.valid_SEQ_structure);
    hm_trimmed = [seqprofile(CB1(selector), 'Alphabet', 'NT', 'Gaps', 'all'), ...
                  seqprofile(CB2(selector), 'Alphabet', 'NT', 'Gaps', 'all'), ...
                  seqprofile(UMI(selector), 'Alphabet', 'NT', 'Gaps', 'all')];
    selector = FQ.masks.valid_lines;
    hm_valid = [seqprofile(CB1(selector), 'Alphabet', 'NT', 'Gaps', 'all'), ...
                seqprofile(CB2(selector), 'Alphabet', 'NT', 'Gaps', 'all'), ...
                seqprofile(UMI(selector), 'Alphabet', 'NT', 'Gaps', 'all')];
            
    [CB1, CB2] = cellfun(@(x) deal(x(1:end-8), x(end-7:end)), vertcat(summary.allele_colony{:}), 'un', false);
    assert(size(CB1,1)==summary.N.called_tags);
    CB1 = pad(CB1, maxlen, 'left', '-');
    dedup = unique([dat.CB(dat.allele>0) dat.UMI(dat.allele>0)], 'rows');
    UMI = FQ.UMI(dedup(:,2));

    hm_cell  = [seqprofile(CB1, 'Alphabet', 'NT', 'Gaps', 'all'), ...
                seqprofile(CB2, 'Alphabet', 'NT', 'Gaps', 'all'), ...
                seqprofile(UMI, 'Alphabet', 'NT', 'Gaps', 'all')];
            
    hm       = cumsum([hm_all; hm_trimmed; hm_valid; hm_cell], 1, 'reverse');
    
    sp = subplot(Nr,Nc,which_sp);
    hold on;
    for i = 1:4        
        bar(hm((i-1)*5+1,:), 'blue',   'BarWidth', 0.99);
        bar(hm((i-1)*5+2,:), 'green',  'BarWidth', 0.99);
        bar(hm((i-1)*5+3,:), 'yellow', 'BarWidth', 0.99);
        bar(hm((i-1)*5+4,:), 'red',    'BarWidth', 0.99);
        bar(hm((i-1)*5+5,:), 'white',  'BarWidth', 0.99);
    end
    width = size(hm,2);
    UMI_width = length(UMI{1});
    patch([width-UMI_width-7-0.5 width-UMI_width-7-0.5 width-UMI_width+0.5 width-UMI_width+0.5], ...
          [0 4 4 0], 'white', 'FaceAlpha', 0.5);
    hold off;
    
    set(gca, 'Xtick', [], 'Ytick', []);
    axis tight; box on;
        
end

function sp = UMI_composition_by_allele(Nr, Nc, which_sp, summary)

    hm = cellfun(@(x) sum(seqprofile(x, 'Alphabet', 'NT'),2), summary.allele_colony, 'un', false);    
    hm = horzcat(hm{:});
    hm = hm./sum(hm,1);
    hm = cumsum(hm, 1, 'reverse');
    N = size(hm, 2);

    sp = subplot(Nr,Nc,which_sp);
    
    if (N>0)
        bar(1:N, hm(1,:), 'blue', 'BarWidth', 0.99);   hold on;
        bar(1:N, hm(2,:), 'green', 'BarWidth', 0.99);
        bar(1:N, hm(3,:), 'yellow', 'BarWidth', 0.99);
        bar(1:N, hm(4,:), 'red', 'BarWidth', 0.99);    hold off;
    end
    set(gca,'YAxisLocation', 'right', 'Ytick', [0.125:0.25:1.0], 'Yticklabels', {'T'; 'G'; 'C'; 'A'}, 'Xtick', []);
    axis tight; box on;

end

function sp = CB_UMI_composition_by_allele(Nr, Nc, which_sp, summary, FQ, dat)
    
    dedup = unique([dat.allele(dat.allele>0) dat.CB(dat.allele>0) dat.UMI(dat.allele>0)], 'rows');
    
    CB1 = cellfun(@(y) cellfun(@(x) x(1:end-8), y, 'un', false), summary.allele_colony, 'un', false);
    maxlen = max(cellfun(@(x) max(cellfun(@length, x)), CB1));
    CB1 = cellfun(@(x) pad(x, maxlen, 'left', '-'), CB1, 'un', false);
    CB2 = cellfun(@(y) cellfun(@(x) x(end-7:end), y, 'un', false), summary.allele_colony, 'un', false);
    UMI = arrayfun(@(i) FQ.UMI(dedup(dedup(:,1)==i,3)), [1:size(summary.allele_colony,1)]', 'un', false);

    hm_cb1 = cellfun(@(x) sum(seqprofile(x, 'Alphabet', 'NT', 'Gaps', 'all'),2), CB1, 'un', false);    
    hm_cb1 = horzcat(hm_cb1{:});
    hm_cb1 = hm_cb1./sum(hm_cb1,1);
    hm_cb2 = cellfun(@(x) sum(seqprofile(x, 'Alphabet', 'NT', 'Gaps', 'all'),2), CB2, 'un', false);    
    hm_cb2 = horzcat(hm_cb2{:});
    hm_cb2 = hm_cb2./sum(hm_cb2,1);
    
    hm_umi = cellfun(@(x) sum(seqprofile(x, 'Alphabet', 'NT', 'Gaps', 'all'),2), UMI, 'un', false);
    hm_umi = horzcat(hm_umi{:});
    hm_umi = hm_umi./sum(hm_umi,1);

    hm = cumsum([hm_umi; hm_cb2; hm_cb1], 1, 'reverse');
    N = size(hm, 2);

    sp = subplot(Nr,Nc,which_sp);

    hold on;
    for i = 1:3
        bar(1:N, hm((i-1)*5+1,:), 'blue',   'BarWidth', 0.99);
        bar(1:N, hm((i-1)*5+2,:), 'green',  'BarWidth', 0.99);
        bar(1:N, hm((i-1)*5+3,:), 'yellow', 'BarWidth', 0.99);
        bar(1:N, hm((i-1)*5+4,:), 'red',    'BarWidth', 0.99);   
        bar(1:N, hm((i-1)*5+5,:), 'white',  'BarWidth', 0.99);   
    end
    hold off;
    
    set(gca,'YAxisLocation', 'right', 'Ytick', [0.1:0.2:3.0], ...
        'Yticklabels', {'-1'; 'T1'; 'G1'; 'C1'; 'A1'; '-2'; 'T2'; 'G2'; 'C2'; 'A2'; '-U'; 'TU'; 'GU'; 'CU'; 'AU'}, 'Xtick', []);
    yrule = get(gca, 'YAxis');
    yrule.FontSize = 6;        
    axis tight; box on;
    
end

function sp = reads_by_tag(Nr, Nc, which_sp, dat, cutoff, summary, type)

    template_index = find(strcmp(cellfun(@(x) degap(x.get_seq), summary.alleles, 'un', false), ...
                                 CARLIN_def.getInstance.seq.CARLIN));

    y_reads = accumarray(dat.CB,1);
    y_reads = y_reads(1:cutoff);

    [x_reads_per_CB, ~, y_reads_per_CB] = unique([dat.CB dat.CB_raw], 'rows');
    y_reads_per_CB = accumarray(y_reads_per_CB,1);
    x_reads_per_CB = x_reads_per_CB(:,1);
    y_reads_per_CB(x_reads_per_CB>cutoff)=[];
    x_reads_per_CB(x_reads_per_CB>cutoff)=[];
    
    [x_reads_called, ~, y_reads_called] = unique(dat.CB(dat.allele~=0));
    y_reads_called = accumarray(y_reads_called,1);
    y_reads_called(x_reads_called > cutoff) = [];
    x_reads_called(x_reads_called > cutoff) = [];
    if (~isempty(template_index))
        template_mask = arrayfun(@(i) any(dat.allele(dat.CB==i)==template_index), x_reads_called);
    else
        template_mask = false(size(x_reads_called));
    end
    
    [x_reads_sanitized, ~, y_reads_sanitized] = unique([dat.CB dat.sanitized dat.allele], 'rows');
    y_reads_sanitized = accumarray(y_reads_sanitized,1);
    y_reads_sanitized(x_reads_sanitized(:,1)>cutoff)  = [];
    x_reads_sanitized(x_reads_sanitized(:,1)>cutoff,:)= [];
    discard_mask      = x_reads_sanitized(:,3)==0;
    x_reads_sanitized = x_reads_sanitized(:,1);
    
    sp = subplot(Nr,Nc,which_sp);
    
    plot(log10(y_reads), 'Color', 'black', 'LineWidth', 2.0); hold on;
    legend_entries = {'Valid Lines'};
    if (any(template_mask))
        y_reads_template = zeros(size(y_reads_called));
        y_reads_template(template_mask) = y_reads_called(template_mask);
        bar(x_reads_called, log10(y_reads_template), ...
            'EdgeColor', 'none', 'FaceColor', 'magenta', 'FaceAlpha', 0.2, 'BarWidth', 0.99);
        legend_entries = [legend_entries; cellstr('Template Allele')];
    end
    if (any(~template_mask))
        y_reads_edited = zeros(size(y_reads_called));
        y_reads_edited(~template_mask) = y_reads_called(~template_mask);
        bar(x_reads_called, log10(y_reads_edited), ...
            'EdgeColor', 'none', 'FaceColor', 'cyan', 'FaceAlpha', 0.4, 'BarWidth', 0.99);    
        legend_entries = [legend_entries; cellstr('Edited Allele')];
    end
    scatter(x_reads_sanitized( discard_mask), log10(y_reads_sanitized( discard_mask)), [], 'b', 'Marker', '.');
    scatter(x_reads_sanitized(~discard_mask), log10(y_reads_sanitized(~discard_mask)), 10, 'g', 'filled', 'Marker', 'o');
    scatter(x_reads_per_CB, log10(y_reads_per_CB), [], 'r', 'Marker', '.');
    hold off;
    
    legend([legend_entries; {'Sanitized - Discarded'; 'Sanitized - Merged in Allele'; sprintf('Per %s', type)}]);
    
    ylabel('log 10 (Reads)');
    
    xlim([0, cutoff]); axis tight; box on;
    set(gca, 'Xtick', []); 
    
end

function sp = reads_by_allele(Nr, Nc, which_sp, dat, summary)
                            
    assert(summary.N.called_tags == length(unique(dat.CB(dat.allele>0))));
    
    template_mask = strcmp(cellfun(@(x) degap(x.get_seq), summary.alleles, 'un', false), ...
                            CARLIN_def.getInstance.seq.CARLIN);
    
    y_reads = accumarray(dat.allele_eventual(dat.allele_eventual>0),1);    
    assert(sum(y_reads) == summary.reads.called_tags_total);
    assert(sum(y_reads(~template_mask)) == summary.reads.eventful_tags_total);
    y_reads_called = accumarray(dat.allele(dat.allele>0),1);
    assert(sum(y_reads_called) == summary.reads.called_tags_allele);
    assert(sum(y_reads_called(~template_mask)) == summary.reads.eventful_tags_allele);
    
    [x_reads_per_CB, ~, y_reads_per_CB] = unique([dat.allele(dat.allele>0) dat.CB(dat.allele>0)], 'rows');
    y_reads_per_CB = accumarray(y_reads_per_CB,1);
    x_reads_per_CB = x_reads_per_CB(:,1);    
        
    [x_reads_sanitized, ~, y_reads_sanitized] = unique(...
        [dat.allele_eventual(dat.allele_eventual>0), dat.allele(dat.allele_eventual>0), dat.sanitized(dat.allele_eventual>0)], 'rows');
    y_reads_sanitized = accumarray(y_reads_sanitized,1);
    discard_mask      = x_reads_sanitized(:,2)==0;
    x_reads_sanitized = x_reads_sanitized(:,1);
    
    sp = subplot(Nr,Nc,which_sp);
    set(gca,'YDir','normal', 'Yticklabels', []); yyaxis right;
    
    plot(log10(y_reads), 'Color', 'black', 'LineWidth', 2.0); hold on;
    if (any(template_mask))
        bar(find(template_mask), log10(y_reads_called (template_mask)), ...
            'EdgeColor', 'none', 'FaceColor', 'magenta', 'FaceAlpha', 0.2, 'BarWidth', 0.99);
    end
    
    if (any(~template_mask))        
        bar(find(~template_mask), log10(y_reads_called(~template_mask)), ...
            'EdgeColor', 'none', 'FaceColor', 'cyan', 'FaceAlpha', 0.4, 'BarWidth', 0.99);    
    end
    scatter(x_reads_sanitized( discard_mask), log10(y_reads_sanitized( discard_mask)), [], 'b', 'Marker', '.');
    scatter(x_reads_sanitized(~discard_mask), log10(y_reads_sanitized(~discard_mask)), 10, 'g', 'filled', 'Marker', 'o');
    scatter(x_reads_per_CB, log10(y_reads_per_CB), [], 'r', 'Marker', '.');
    hold off;
        
    axis tight; box on;
    set(gca, 'Xtick', []); 
 
end
                 
function sp = distance_by_tag(Nr, Nc, which_sp, dat, cutoff, FQ, aligned, summary)
                     
    y_sanitized = unique([dat.CB dat.sequence dat.sanitized], 'rows');
    y_sanitized(y_sanitized(:,1)>cutoff,:) = [];
    
    L_SEQ_valid   = cellfun(@length, FQ.SEQ_valid);
    L_SEQ_aligned = cellfun(@(x) length(degap(x.get_seq())), aligned.aligned_SEQ);
    
    bp_diff = abs(L_SEQ_valid(y_sanitized(:,2))-L_SEQ_aligned(y_sanitized(:,3)));
    
    [y_sanitized_unique, ~, y_sanitized_ind] = unique([y_sanitized(bp_diff==0,2) y_sanitized(bp_diff==0,3)], 'rows');
    temp_bp_diff_0 = cellfun(@(x,y) sum(int2nt(x)~=y), FQ.SEQ_valid(y_sanitized_unique(:,1)), ...
                             cellfun(@(z) degap(z.get_seq()), aligned.aligned_SEQ(y_sanitized_unique(:,2)), 'un', false));
    bp_diff(bp_diff==0) = temp_bp_diff_0(y_sanitized_ind);
    
%     
%     bp_diff = cellfun(@(x,y) abs(length(x)-length(y)), FQ.SEQ_valid(y_sanitized(:,2)), ...
%              cellfun(@(z) degap(z.get_seq()), aligned.aligned_SEQ(y_sanitized(:,3)), 'un', false));
% 
%     bp_diff(bp_diff==0) = cellfun(@(x,y) sum(int2nt(x)~=y), FQ.SEQ_valid(y_sanitized(bp_diff==0,2)), ...
%                           cellfun(@(z) degap(z.get_seq()), aligned.aligned_SEQ(y_sanitized(bp_diff==0,3)), 'un', false));

    x_sanitized = unique([dat.CB dat.sanitized dat.allele], 'rows');
    x_sanitized(x_sanitized(:,1)>cutoff,:) = [];
    discard_mask = x_sanitized(:,3)==0;
    
%     L_sanitized = cellfun(@(x) length(degap(x.get_seq)), aligned.aligned_SEQ(x_sanitized(:,2)));
    L_sanitized     = L_SEQ_aligned(x_sanitized(:,2));
    L_sanitized_neg = zeros(size(L_sanitized));
    L_sanitized_pos = zeros(size(L_sanitized));
    
    % Don't bother computing edit distances for sanitized sequences that
    % end up being discarded being this really clutters the plot
    [L_sanitized_neg(~discard_mask), L_sanitized_pos(~discard_mask)] = ...
        arrayfun(@(i,j) deal(min(bp_diff(y_sanitized(:,1)==i & y_sanitized(:,3)==j)), ...
                             max(bp_diff(y_sanitized(:,1)==i & y_sanitized(:,3)==j))), ...
                        x_sanitized(~discard_mask,1), x_sanitized(~discard_mask,2));
    
    y_allele = unique([dat.CB dat.sanitized dat.allele], 'rows');
    y_allele = y_allele(y_allele(:,3)>0,:);
    y_allele(y_allele(:,1)>cutoff,:) = [];
    
    bp_diff = cellfun(@(x,y) abs(length(x)-length(y)), ...
              cellfun(@(z) degap(z.get_seq()), aligned.aligned_SEQ(y_allele(:,2)), 'un', false), ...
              cellfun(@(z) degap(z.get_seq()), summary.alleles(y_allele(:,3)), 'un', false));
    
    bp_diff(bp_diff==0) = cellfun(@(x,y) sum(x~=y), ...
                          cellfun(@(z) degap(z.get_seq()), aligned.aligned_SEQ(y_allele(bp_diff==0,2)), 'un', false), ...
                          cellfun(@(z) degap(z.get_seq()), summary.alleles(y_allele(bp_diff==0,3)), 'un', false));
                  
    x_allele = unique([dat.CB(dat.allele>0), dat.allele(dat.allele>0)], 'rows');
    x_allele(x_allele(:,1)>cutoff,:) = [];
    
    L_allele = cellfun(@(x) length(degap(x.get_seq)), summary.alleles(x_allele(:,2)));
    [L_allele_neg, L_allele_pos] = arrayfun(@(i) deal(min(bp_diff(y_allele(:,1)==i)), ...
                                                      max(bp_diff(y_allele(:,1)==i))), x_allele(:,1));
        
    sp = subplot(Nr,Nc,which_sp);
    
    s = scatter(x_sanitized( discard_mask,1), L_sanitized( discard_mask), [], 'b', 'Marker', '.'); hold on;
    s.MarkerFaceAlpha = 0.3;
    s.MarkerEdgeAlpha = 0.3;
    scatter(x_sanitized(~discard_mask,1), L_sanitized(~discard_mask), [], 'g', 'Marker', '.'); 
    scatter(x_allele(:,1),                L_allele,                   [], 'r', 'Marker', '.'); 
%     errorbar(x_sanitized(discard_mask,1), L_sanitized(discard_mask), ...
%              L_sanitized_neg(discard_mask), L_sanitized_pos(discard_mask), ...
%              'LineStyle', 'None', 'color', 'b');
    errorbar(x_sanitized(~discard_mask,1), L_sanitized(~discard_mask), ...
             L_sanitized_neg(~discard_mask), L_sanitized_pos(~discard_mask), ...
             'LineStyle', 'None', 'color', 'g'); 
    errorbar(x_allele(:,1), L_allele, L_allele_neg, L_allele_pos, 'LineStyle', 'None', 'color', 'r'); hold off;
    
    legend({'Sanitized - Discarded'; 'Sanitized - Merged in Allele'; 'Allele'});
    ylabel('Length +/- Edit Distance (bp)');
    
    xlim([0, cutoff]); axis tight; box on;
    set(gca, 'Xtick', []); 
end

function [sp, suspect_alleles] = distance_by_allele(Nr, Nc, which_sp, dat, FQ, aligned, summary, type)

    if (strcmp(type, 'UMI'))
        
        sp = subplot(Nr,Nc,which_sp);
        
        if (~isempty(summary.allele_colony))
            L = length(summary.allele_colony{1}{1});
            hm = cellfun(@(x) histcounts(seqpdist(x, 'Method', 'p-distance', 'Alphabet', 'NT', 'UseParallel', true)*L,...
                                         -0.5:1.0:L+0.5)', summary.allele_colony, 'un', false);
            hm = horzcat(hm{:});
            hm = hm./sum(hm,1);
            hm(isnan(hm)) = 0;

            suspect_alleles = find(sum(hm(1:2,:),1) >= 0.1);            
        else
            hm = 1;
            suspect_alleles = [];
        end
        imagesc(1-hm); colormap gray; caxis([0, 1]); grid on;
        set(gca,'YDir','normal', 'Yticklabels', [], 'Xtick', []);
    else
        
        sp = subplot(Nr,Nc,which_sp);
        
        if (~isempty(summary.allele_colony))
            dedup = unique([dat.allele(dat.allele>0) dat.CB(dat.allele>0) dat.UMI(dat.allele>0)], 'rows');

            CB1 = cellfun(@(y) cellfun(@(x) x(1:end-8), y, 'un', false), summary.allele_colony, 'un', false);
            maxlen = max(cellfun(@(x) max(cellfun(@length, x)), CB1));
            CB1 = cellfun(@(x) pad(x, maxlen, 'left', '-'), CB1, 'un', false);
            CB2 = cellfun(@(y) cellfun(@(x) x(end-7:end), y, 'un', false), summary.allele_colony, 'un', false);
            UMI = arrayfun(@(i) FQ.UMI(dedup(dedup(:,1)==i,3)), [1:size(summary.allele_colony,1)]', 'un', false);

            hm1 = cellfun(@(x) histcounts(seqpdist(x, 'Method', 'p-distance', 'Alphabet', 'NT')*maxlen,-0.5:1.0:maxlen+0.5)', ...
                              CB1, 'un', false);
            hm1 = horzcat(hm1{:});
            hm1 = hm1./sum(hm1,1);
            hm1(isnan(hm1)) = 0;

            hm2 = cellfun(@(x) histcounts(seqpdist(x, 'Method', 'p-distance', 'Alphabet', 'NT')*8,-0.5:1.0:8+0.5)', ...
                              CB2, 'un', false);
            hm2 = horzcat(hm2{:});
            hm2 = hm2./sum(hm2,1);
            hm2(isnan(hm2)) = 0;

            L_UMI = length(UMI{1}{1});
            hm3 = cellfun(@(x) histcounts(seqpdist(x, 'Method', 'p-distance', 'Alphabet', 'NT')*L_UMI,-0.5:1.0:L_UMI+0.5)', ...
                              UMI, 'un', false);
            hm3 = horzcat(hm3{:});
            hm3 = hm3./sum(hm3,1);
            hm3(isnan(hm3)) = 0;

            suspect_alleles = find(sum(hm1(1:2,:),1)>=0.1 | sum(hm2(1:2,:),1)>=0.1 | sum(hm3(1:2,:))>=0.1);

            hm = [hm1; hm2; hm3];

            imagesc(1-hm); colormap gray; caxis([0, 1]); grid on; hold on;
            patch([0.5, 0.5, size(hm,2)+0.5 size(hm,2)+0.5], ...
                  [maxlen+2-0.5 maxlen+2+8+0.5 maxlen+2+8+0.5 maxlen+2-0.5], ...
                  'yellow', 'FaceAlpha', 0.2, 'EdgeColor', 'yellow');
            hold off;
            
        else
            hm = 1;
            suspect_alleles = [];
            imagesc(1-hm); colormap gray; caxis([0, 1]); grid on; hold on;
        end
        set(gca,'YDir','normal', 'Yticklabels', [], 'Xtick', []);
     
    end

    y_sanitized = unique([dat.allele_eventual dat.sequence dat.sanitized], 'rows');
    y_sanitized(y_sanitized(:,1)==0,:) = [];
     
    L_SEQ_valid   = cellfun(@length, FQ.SEQ_valid);
    L_SEQ_aligned = cellfun(@(x) length(degap(x.get_seq())), aligned.aligned_SEQ);
    
    bp_diff = abs(L_SEQ_valid(y_sanitized(:,2))-L_SEQ_aligned(y_sanitized(:,3)));

    [y_sanitized_unique, ~, y_sanitized_ind] = unique([y_sanitized(bp_diff==0,2) y_sanitized(bp_diff==0,3)], 'rows');
    temp_bp_diff_0 = cellfun(@(x,y) sum(int2nt(x)~=y), FQ.SEQ_valid(y_sanitized_unique(:,1)), ...
                             cellfun(@(z) degap(z.get_seq()), aligned.aligned_SEQ(y_sanitized_unique(:,2)), 'un', false));
    bp_diff(bp_diff==0) = temp_bp_diff_0(y_sanitized_ind);
% 
%     
%     bp_diff(bp_diff==0) = cellfun(@(x,y) sum(int2nt(x)~=y), FQ.SEQ_valid(y_sanitized(bp_diff==0,2)), ...
%                           cellfun(@(z) degap(z.get_seq()), aligned.aligned_SEQ(y_sanitized(bp_diff==0,3)), 'un', false));

    x_sanitized = unique([dat.allele_eventual dat.allele dat.sanitized], 'rows');
    x_sanitized(x_sanitized(:,1)==0,:) = [];
    discard_mask = x_sanitized(:,2)==0;
    
%     L_sanitized = cellfun(@(x) length(degap(x.get_seq)), aligned.aligned_SEQ(x_sanitized(:,3)));
    L_sanitized     = L_SEQ_aligned(x_sanitized(:,3));
    L_sanitized_neg = zeros(size(L_sanitized));
    L_sanitized_pos = zeros(size(L_sanitized));
    
    % Don't bother computing edit distances for sanitized sequences that
    % end up being discarded because this really clutters the plot    
    [L_sanitized_neg(~discard_mask), L_sanitized_pos(~discard_mask)] = ...
        arrayfun(@(i,j) deal(min(bp_diff(y_sanitized(:,1)==i & y_sanitized(:,3)==j)), ...
                             max(bp_diff(y_sanitized(:,1)==i & y_sanitized(:,3)==j))), ...
                        x_sanitized(~discard_mask,1), x_sanitized(~discard_mask,3));
    
    x_sanitized = x_sanitized(:,1);
    
    y_allele = unique([dat.allele dat.sanitized], 'rows');
    y_allele(y_allele(:,1)==0,:) = [];
    
    bp_diff = cellfun(@(x,y) abs(length(x)-length(y)), ...
                      cellfun(@(z) degap(z.get_seq()), summary.alleles(y_allele(:,1)), 'un', false), ...
                      cellfun(@(z) degap(z.get_seq()), aligned.aligned_SEQ(y_allele(:,2)), 'un', false));
    
    bp_diff(bp_diff==0) = cellfun(@(x,y) sum(x~=y), ...
                          cellfun(@(z) degap(z.get_seq()), summary.alleles(y_allele(bp_diff==0,1)), 'un', false), ...
                          cellfun(@(z) degap(z.get_seq()), aligned.aligned_SEQ(y_allele(bp_diff==0,2)), 'un', false));
                      
                  
    x_allele = unique(dat.allele(dat.allele>0));
    L_allele = cellfun(@(x) length(degap(x.get_seq)), summary.alleles(x_allele));
    [L_allele_neg, L_allele_pos] = arrayfun(@(i) deal(min(bp_diff(y_allele(:,1)==i)), ...
                                                      max(bp_diff(y_allele(:,1)==i))), x_allele);
    
    set(gca,'YDir','normal', 'Yticklabels', []); yyaxis right;

    s = scatter(x_sanitized( discard_mask), L_sanitized( discard_mask), [], 'b', 'Marker', '.'); hold on;
    s.MarkerFaceAlpha = 0.3;
    s.MarkerEdgeAlpha = 0.3;
    scatter(x_sanitized(~discard_mask), L_sanitized(~discard_mask), [], 'g', 'Marker', '.'); 
    scatter(x_allele, L_allele, [], 'r', 'Marker', '.'); 

%     errorbar(x_sanitized(discard_mask), L_sanitized(discard_mask), ...
%              L_sanitized_neg(discard_mask), L_sanitized_pos(discard_mask), 'LineStyle', 'None', 'Color', 'blue');
    errorbar(x_sanitized(~discard_mask), L_sanitized(~discard_mask), ...
             L_sanitized_neg(~discard_mask), L_sanitized_pos(~discard_mask), 'LineStyle', 'None', 'Color', 'green');         
    errorbar(x_allele, L_allele, L_allele_neg, L_allele_pos, 'LineStyle', 'None', 'Color', 'red'); 
    hold off;
    
    axis tight; box on; 
    set(gca, 'Xtick', []);
end

function sp = variety_by_tag(Nr, Nc, which_sp, dat, cutoff, type)

    num_CBs        = unique([dat.CB dat.CB_raw  ], 'rows');    
    num_sequences   = unique([dat.CB dat.sequence ], 'rows');        
    num_sanitized   = unique([dat.CB dat.sanitized], 'rows');    
    num_events      = unique([dat.CB dat.event    ], 'rows');    

    num_CBs        = accumarray(num_CBs     (:,1), 1);
    num_sequences   = accumarray(num_sequences(:,1), 1);
    num_sanitized   = accumarray(num_sanitized(:,1), 1);
    num_events      = accumarray(num_events   (:,1), 1);
    
    num_CBs        = num_CBs(1:cutoff);
    num_sequences   = num_sequences(1:cutoff);
    num_sanitized   = num_sanitized(1:cutoff);
    num_events      = num_events(1:cutoff);
    
    if (strcmp(type, 'CB'))
        num_UMIs     = unique([dat.CB dat.UMI],     'rows');
        num_UMIs_raw = unique([dat.CB dat.UMI_raw], 'rows');
        num_UMIs     = accumarray(num_UMIs    (:,1), 1);
        num_UMIs_raw = accumarray(num_UMIs_raw(:,1), 1);
        num_UMIs     = num_UMIs(1:cutoff);
        num_UMIs_raw = num_UMIs_raw(1:cutoff);
    end
    
    sp = subplot(Nr,Nc,which_sp);

    bar(log10(num_sequences), 'FaceColor', [0.5 0.5 1.0], 'BarWidth', 0.99); hold on;
    bar(log10(num_sanitized), 'FaceColor', [0.5 1.0 0.5], 'BarWidth', 0.99);
    bar(log10(num_events),    'FaceColor', [1.0 0.5 0.5], 'BarWidth', 0.99);
    plot(log10(num_CBs),      'Marker', '.', 'LineStyle', 'None', 'Color', 'black'); 
    if (strcmp(type, 'CB'))
        plot(log10(num_UMIs_raw),   'Marker', '.', 'LineStyle', 'None', 'Color', 'magenta');
        plot(log10(num_UMIs),       'Marker', '.', 'LineStyle', 'None', 'Color', 'cyan');
    end
    hold off;    
    if (strcmp(type, 'UMI'))
        legend({'Sequence'; 'Sanitized'; 'Structure'; 'Noisy UMIs'});
    elseif (strcmp(type, 'CB'))
        legend({'Sequence'; 'Sanitized'; 'Structure'; 'Noisy CBs'; 'Noisy UMIs'; 'UMIs'});
    end
    ylabel('log10');
    
    xlim([0, cutoff]); axis tight; box on;
    xlabel(sprintf('%s (ranked by valid reads)', type));
    
end

function sp = variety_by_allele(Nr, Nc, which_sp, dat, type)

    num_CBs               = unique([dat.allele dat.CB               ], 'rows');
    num_sequences_duped   = unique([dat.allele dat.CB dat.sequence  ], 'rows');
    num_sequences_deduped = unique([dat.allele dat.sequence         ], 'rows');
    num_sanitized_duped   = unique([dat.allele dat.CB dat.sanitized ], 'rows');
    num_sanitized_deduped = unique([dat.allele dat.sanitized        ], 'rows');
        
    num_CBs              (num_CBs              (:,1)==0,:) = [];
    num_sequences_duped  (num_sequences_duped  (:,1)==0,:) = [];
    num_sequences_deduped(num_sequences_deduped(:,1)==0,:) = [];
    num_sanitized_duped  (num_sanitized_duped  (:,1)==0,:) = [];
    num_sanitized_deduped(num_sanitized_deduped(:,1)==0,:) = [];

    num_CBs               = accumarray(num_CBs              (:,1), 1);
    num_sequences_duped   = accumarray(num_sequences_duped  (:,1), 1);
    num_sequences_deduped = accumarray(num_sequences_deduped(:,1), 1);
    num_sanitized_duped   = accumarray(num_sanitized_duped  (:,1), 1);
    num_sanitized_deduped = accumarray(num_sanitized_deduped(:,1), 1);
    
    if (strcmp(type, 'CB'))
        num_UMIs_duped = unique([dat.allele dat.CB dat.UMI], 'rows');
        num_UMIs_duped(num_UMIs_duped(:,1)==0,:) = [];
        num_UMIs_duped = accumarray(num_UMIs_duped(:,1),1);
        num_UMIs_deduped = unique([dat.allele dat.UMI], 'rows');
        num_UMIs_deduped(num_UMIs_deduped(:,1)==0,:) = [];
        num_UMIs_deduped = accumarray(num_UMIs_deduped(:,1),1);
    end 
    
    sp = subplot(Nr,Nc,which_sp);
    set(gca,'Ytick', [], 'Yticklabels', []); yyaxis right;

    bar(log10(num_sequences_duped),   'FaceColor', 'b', 'BarWidth', 0.99); hold on;
    bar(log10(num_sequences_deduped), 'FaceColor', 'c', 'BarWidth', 0.99); 
    bar(log10(num_sanitized_duped),   'FaceColor', 'g', 'BarWidth', 0.99);
    bar(log10(num_sanitized_deduped), 'FaceColor', 'y', 'BarWidth', 0.99);
    if (strcmp(type, 'CB'))
        plot(log10(num_UMIs_duped), 'Color', 'red', 'LineStyle', '-', 'LineWidth', 2.0, 'Marker', 'o');
        plot(log10(num_UMIs_deduped), 'Color', [0.5, 0 0], 'LineStyle', '-', 'LineWidth', 2.0, 'Marker', 'd');
    end
    plot(log10(num_CBs), 'Color', 'black', 'LineStyle', '-', 'LineWidth', 2.0, 'Marker', 's');
    hold off;
    
    if (strcmp(type, 'UMI'))
        legend({'Sequences (UMI dup)'; 'Sequences (UMI dedup)'; ...
                'Sanitized (UMI dup)'; 'Sanitized (UMI dedup)'; 'UMIs'});
        xlabel('Allele (Ranked by UMI Freq)');
    else
        legend({'Sequences (CB dup)'; 'Sequences (CB dedup)'; ...
                'Sanitized (CB dup)'; 'Sanitized (CB dedup)'; 'UMIs (CB dup)'; 'UMIs (CB dedup)'; 'CBs'});
        xlabel('Allele (Ranked by CB Freq)');
    end
    axis tight; box on;
    
end

function position_subplots(N1, N2, sp)
    
    tight_margin = 0.01;
    bottom_label = 0.08;
    top_label    = 0.06;
    left_label = 0.04;    
    right_label = 0.02;
    rel_height = (1-bottom_label-3*tight_margin)/4;

    min_left_panel = 0.4;
    min_right_panel = 0.2;    
    divvy_up = N1/(N1+N2);
    width_leftover = 1-tight_margin-left_label-right_label-min_left_panel-min_right_panel;    
    rel_left = min_left_panel+divvy_up*width_leftover;
    rel_right = min_right_panel+(1-divvy_up)*width_leftover;
   
    top_left_label = 0.07;
    rel_left_panel1 = 0.1;    
    rel_left_panel2 = 0.1;
    rel_left_panel4 = 0.1;    
    rel_left_top_widespace = 0.02;
    
    rel_left_panel3   = (left_label+rel_left)-...
                            (top_left_label+rel_left_panel1+rel_left_panel2+rel_left_panel4+ ...
                            2*rel_left_top_widespace+tight_margin);
                        
    set(sp(4,1), 'Position',[left_label, ...
                             bottom_label+0*rel_height, ...
                             rel_left, ...
                             rel_height]);
                         
    set(sp(4,2), 'Position',[left_label+rel_left+tight_margin, ...
                             bottom_label+0*rel_height, ...
                             rel_right, ...
                             rel_height]);
    
    set(sp(3,1), 'Position',[left_label, ...
                             bottom_label+1*rel_height+1*tight_margin, ...
                             rel_left, ...
                             rel_height]);
                         
    set(sp(3,2), 'Position',[left_label+rel_left+tight_margin, ...
                             bottom_label+1*rel_height+1*tight_margin, ...
                             rel_right, ...
                             rel_height]);
    
    set(sp(2,1), 'Position',[left_label, ...
                             bottom_label+2*rel_height+2*tight_margin, ...
                             rel_left, ...
                             rel_height]);
                         
    set(sp(2,2), 'Position',[left_label+rel_left+tight_margin, ...
                             bottom_label+2*rel_height+2*tight_margin, ...
                             rel_right, ...
                             rel_height]);    

    set(sp(1,1), 'Position',[top_left_label, ...
                             bottom_label+3*rel_height+3*tight_margin, ...
                             rel_left_panel1, ...
                             rel_height-top_label]);
                         
    set(sp(1,2), 'Position',[top_left_label+rel_left_panel1+rel_left_top_widespace, ...
                             bottom_label+3*rel_height+3*tight_margin, ...
                             rel_left_panel1, ...
                             rel_height-top_label]);
                         
    set(sp(1,3), 'Position',[top_left_label+rel_left_panel1+rel_left_panel2+rel_left_top_widespace+tight_margin, ...
                             bottom_label+3*rel_height+3*tight_margin, ...
                             rel_left_panel3, ...
                             rel_height-top_label]);
                         
    set(sp(1,4), 'Position',[top_left_label+rel_left_panel1+rel_left_panel2+rel_left_panel3+tight_margin+2*rel_left_top_widespace, ...
                             bottom_label+3*rel_height+3*tight_margin, ...
                             rel_left_panel4, ...
                             rel_height-top_label]);
                         
    set(sp(1,5), 'Position',[left_label+rel_left+tight_margin, ...
                             bottom_label+3*rel_height+3*tight_margin, ...
                             rel_right, ...
                             rel_height]);     
end