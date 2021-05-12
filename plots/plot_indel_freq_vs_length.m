function [sp, ins_freq, del_freq] = plot_indel_freq_vs_length(summary, allele_breakdown_by_sample, Nr, Nc, which_sp)

    assert(isa(summary, 'ExperimentSummary'));  
    
    if (nargin == 1)
        allele_breakdown_by_sample = summary.allele_freqs;
    end
    
    N_samples = size(allele_breakdown_by_sample,2);
    assert(N_samples <= 3, 'Extend alpha values to use more samples');

    if (nargin < 3)
        if (N_samples == 3)  
            fig_height = 5;
            fig_width = 8.4;
            top_margin = 0.3;
        elseif (N_samples == 2)            
            fig_height = 3.2;
            fig_width = 8;
            top_margin = 0.3;
        else       
            fig_height= 5.0;
            fig_width  = 5.8;
            top_margin = 0.1;
        end
        left_margin = 0.8;
        bottom_margin = 0.7;        
        right_margin = 0.1;        
        figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
               'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
        sp = subplot(1,1,1);
        set(sp, 'Units', 'centimeters', 'Position', [left_margin bottom_margin fig_width-left_margin-right_margin, fig_height-bottom_margin-top_margin]);        
    else
        sp = subplot(Nr, Nc, which_sp);
    end
    
    mut_events = cellfun(@Mutation.identify_Cas9_events, summary.alleles, 'un', false);
    [num_bp_del, num_bp_ins] = cellfun(@(x) arrayfun(@(i) x(i).num_bps_indel, [1:length(x)]'), mut_events, 'un', false);
    L_max = max([vertcat(num_bp_del{:}); vertcat(num_bp_ins{:})]);
    
    num_bp_del = cellfun(@(x) nonzeros(x), num_bp_del, 'un', false);
    num_bp_ins = cellfun(@(x) nonzeros(x), num_bp_ins, 'un', false);
    N_del = cellfun(@(x) length(x), num_bp_del);
    N_ins = cellfun(@(x) length(x), num_bp_ins);
    sample_mask = arrayfun(@(i) sum(logical(allele_breakdown_by_sample),2)==i, [1:N_samples], 'un', false);
        
    del_freq = arrayfun(@(i) accumarray(vertcat(num_bp_del{sample_mask{i}}), ...
                                        repelem(summary.allele_freqs(sample_mask{i}), N_del(sample_mask{i})),[L_max, 1]), ...
                             [1:N_samples], 'un', false);
    ins_freq = arrayfun(@(i) accumarray(vertcat(num_bp_ins{sample_mask{i}}), ...
                                        repelem(summary.allele_freqs(sample_mask{i}), N_ins(sample_mask{i})),[L_max, 1]), ...
                             [1:N_samples], 'un', false);
                         
    if (N_samples == 1)
        ins_freq = cellfun(@(x) x/summary.N.eventful_tags, ins_freq, 'un', false);
        del_freq = cellfun(@(x) x/summary.N.eventful_tags, del_freq, 'un', false);
    end

    h1 = bar([1:L_max]',  horzcat(ins_freq{:}), 'stacked', 'FaceColor', 'b', 'BarWidth', 1.0); hold on;
    h2 = bar([1:L_max]', -horzcat(del_freq{:}), 'stacked', 'FaceColor', 'r', 'BarWidth', 1.0); hold off;
    axis tight;
    xlim([0.5 L_max+0.5]);
    
    if (N_samples == 1)
        alpha_val = 1.0;
    elseif (N_samples == 2)
        alpha_val = [0.4 1.0];
    elseif (N_samples == 3)
        alpha_val = [0.6 0.2, 1.0];
    end
    
    for i = 1:N_samples
        h1(i).FaceAlpha = alpha_val(i);
        h2(i).FaceAlpha = alpha_val(i);
    end    
    
    set(get(gca, 'YAxis'), 'FontSize', 5);
    set(get(gca, 'XAxis'), 'FontSize', 5);    
    set(gca, 'TickDir', 'out', 'LineWidth', 1.0);
    
    xlabel('Mutation Length (bp)', 'FontSize', 6);
    
    if (N_samples > 1)
        ylabel('# Occurrences in Transcripts', 'FontSize', 6);
    else
        ylabel('Fraction of Edited Cells', 'FontSize', 6);
    end
    
    box off;
    
    if (N_samples > 1)
        ax = gca;
        ax.YAxis.Exponent = floor(log10(max(abs(get(gca, 'Ytick')))));
        ax.YRuler.SecondaryLabel.Visible = 'on';        
    end
    set(gca, 'Yticklabels', strrep(get(gca, 'Yticklabels'), '-', ''));

    
    if (N_samples == 1)    
        lgd = legend({'Insertion'; 'Deletion'}, 'Location', 'Southeast', 'FontSize', 5);
        legend('boxoff');
    else
        lgd = legend([repmat({''}, N_samples, 1); strrep(arrayfun(@(i) sprintf('%d Mice', i), [1:N_samples]', 'un', false), {'1 Mice'}, {'1 Mouse'})], ...
                     'NumColumns', 2, 'Location', 'East', 'FontSize', 5);
        legend('boxoff');
        lgd.Units = 'centimeters';
        leg_width = 4;
        leg_height = 1.0;
        lgd.Position = [fig_width-right_margin-leg_width fig_height-leg_height-2*top_margin leg_width, leg_height];
        label_width = 1.0;
        label_height = 0.3;
        an1 = annotation('textbox', 'EdgeColor', 'none', 'string', 'Insertion', 'Margin', 0, 'FontSize', 6);
        tight_margin =0.1;
        an1.Units = 'centimeters';
        an1.Position = [lgd.Position(1)+3*tight_margin lgd.Position(2)+lgd.Position(4)-tight_margin label_width label_height];
        an2 = annotation('textbox', 'EdgeColor', 'none', 'string', 'Deletion', 'Margin', 0, 'FontSize', 6);
        an2.Units = 'centimeters';
        an2.Position = [lgd.Position(1)+3*tight_margin+lgd.Position(3)*0.4 lgd.Position(2)+lgd.Position(4)-tight_margin label_width label_height];
    end
end