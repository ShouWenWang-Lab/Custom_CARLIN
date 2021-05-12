function plot_summary(summary, outdir)

    assert(isa(summary, 'ExperimentSummary'));

    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    scrsz = get(0,'ScreenSize');
    scrsz = scrsz(3:4);
    
    tight_margin = 0.01;
    rel_freq_width = 0.2;
    rel_edit_width = 1-rel_freq_width-3*tight_margin;    
    rel_event_height = 0.2;
    rel_seq_height = 1-rel_event_height-3*tight_margin;    
   
    % Top Left
    s_tl = subplot(2, 2, 1);
    if (isa(summary, 'ExperimentReport'))
        if (isa(summary, 'BulkExperimentReport'))
            tag = 'UMI';
        else
            tag = 'CB';
        end
        text(0, 0.65, {sprintf('Total reads in FASTQ  : \t%8d', summary.reads.in_fastq);
                       sprintf('Valid reads considered: \t%8d', summary.reads.valid_lines);
                       sprintf('%5d %ss denoised to %5d', summary.N.uncleaned_tags, tag, summary.N.cleaned_tags);
                       sprintf('%5d %8s %ss with\t\t%8d reads', summary.N.common_tags, 'common', tag, summary.reads.common_tags);
                       sprintf('%5d %8s %ss with\t%8d reads', summary.N.called_tags, 'callable', tag, summary.reads.called_tags_total);
                       sprintf('%5d %8s %ss with\t%8d reads', summary.N.eventful_tags, 'eventful', tag, summary.reads.eventful_tags_total);
                       sprintf('Fraction of %ss edited:   \t%5.1f %%', tag, max(double(summary.N.eventful_tags)/summary.N.called_tags*100,0));
                       sprintf('%5d alleles', size(summary.alleles,1))}, 'FontName', 'FixedWidth');
    else
        tag = 'tag';
        refseq = CARLIN_def.getInstance.seq.CARLIN;
        is_template = ismember(cellfun(@(x) degap(x.get_seq), summary.alleles, 'un', false), refseq);
        text(0, 0.85, {sprintf('%5d callable %ss', sum(summary.allele_freqs), tag);
                       sprintf('%5d eventful %ss', sum(summary.allele_freqs(~is_template)), tag);
                       sprintf('%5.1f%% of %ss edited', max(sum(summary.allele_freqs(~is_template))/sum(summary.allele_freqs)*100,0), tag);                       
                       sprintf('%5d alleles', size(summary.alleles,1))}, 'FontName', 'FixedWidth');
    end
    axis off;
    
    % Bottom Left
    s_bl = subplot(2,2,3);    
    if (~isempty(summary.allele_freqs))
        barh(log10(summary.allele_freqs), 'BarWidth', 1); axis tight;        
        yticks(s_bl, unique([find(summary.allele_freqs > 1, 1, 'last'), length(summary.allele_freqs)]));
    end
    set(gca, 'YAxisLocation', 'right', 'XAxisLocation', 'top', 'YDir', 'reverse', 'XDir', 'reverse');
    ylabel('Allele Rank', 'FontSize', 10); 
    xlabel(sprintf('log10 %s', tag), 'FontSize', 10);
      
    % Top Right
    s_tr = subplot(2,2,2);
    if (~isempty(summary.alleles))
        [~, del_event, ins_event] = cellfun(@(x) Mutation.classify_bp_event(x), summary.alleles, 'Un', false);
        del_event = vertcat(del_event{:});
        ins_event = vertcat(ins_event{:});
        del_pct = sum(del_event.*summary.allele_freqs,1)/sum(summary.allele_freqs)*100;
        ins_pct = sum(ins_event.*summary.allele_freqs,1)/sum(summary.allele_freqs)*100;
        ymax = max(max(del_pct), max(ins_pct));
        ylim_max = roundn(ymax, 1);
        if (ylim_max < ymax)
            ylim_max = ylim_max+10;
        end
        ylim_max = max(ylim_max, 10);
        if (isempty(ylim_max))
            ylim_max = 10;
        end
        ref = CARLIN_def.getInstance;
        bp_bounds = ref.bounds.ordered;
        hold on;
        for i = 1:ref.N.motifs
            x = [bp_bounds(i,1)-0.5 bp_bounds(i,1)-0.5 bp_bounds(i,2)+0.5 bp_bounds(i,2)+0.5];
            y = [0 ylim_max, ylim_max, 0];
            if (ismember(i, ref.motifs.prefix))
                face_alpha = 1-ref.alpha.prefix;        
            elseif (ismember(i, ref.motifs.consites))
                face_alpha = 1-ref.alpha.consite;
            elseif (ismember(i, ref.motifs.cutsites))
                face_alpha = 1-ref.alpha.cutsite;
            elseif (ismember(i, ref.motifs.pams))
                face_alpha = 1-ref.alpha.pam;
            elseif (ismember(i, ref.motifs.postfix))
                face_alpha = 1-ref.alpha.postfix;    
            end
            patch(x, y, 'black', 'FaceAlpha', face_alpha*ref.alpha.overlay);
        end

        plot(del_pct, 'red');            
        bar(ins_pct, 'blue'); hold off;
        axis tight;  
        ylim([0 ylim_max]);        
    end
    xticks([]);
    ylabel(sprintf('%% %ss Edited', tag), 'FontSize', 10);
    
    % Bottom Right
    s_br = subplot(2,2,4);                      
    if (~isempty(summary.alleles))
        RGB = get_sequence_coloring(summary.alleles, 'bp');

        imshow(RGB, 'XData', [1 size(RGB,2)]', 'YData', [1 size(RGB,1)]); hold on;
        h = imshow(repmat(CARLIN_def.getInstance.alpha.CARLIN, [size(RGB,1), 1]), 'XData', [1 size(RGB,2)]', 'YData', [1 size(RGB,1)]);  hold off;
        
        set(h, 'AlphaData', CARLIN_def.getInstance.alpha.overlay);
        axis tight;        
        daspect([scrsz(2)*rel_seq_height/length(summary.alleles), scrsz(1)*rel_edit_width/CARLIN_def.getInstance.width.CARLIN, 1]);
    end    
    
    axis on; box on; xticks([]); yticks([]);
    linkaxes([s_bl, s_br], 'y');
    linkaxes([s_br, s_tr], 'x');
    
    set(s_tl, 'Position',[tight_margin 2*tight_margin+rel_seq_height, rel_freq_width, rel_event_height]);
    set(s_tr, 'Position',[2*tight_margin+rel_freq_width+rel_edit_width*0.045, ...
                          2*tight_margin+rel_seq_height, rel_edit_width*0.91, rel_event_height]);
    set(s_bl, 'Position',[tight_margin tight_margin, rel_freq_width, rel_seq_height]);
    set(s_br, 'Position',[2*tight_margin+rel_freq_width+rel_edit_width*0.045, ...
                          tight_margin, rel_edit_width*0.91, rel_seq_height]);
    
    if (nargin > 1)
        if (~exist(outdir, 'dir'))
            mkdir(outdir);
        end
        print(fig, sprintf('%s/Alleles.png', outdir),'-dpng','-r0');
        close;
    end
end