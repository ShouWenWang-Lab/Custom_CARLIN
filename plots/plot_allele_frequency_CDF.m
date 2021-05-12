function sp = plot_allele_frequency_CDF(summary, title_str)

    assert(isa(summary, 'ExperimentSummary'));
    
    fig_width  = 10.8;
    fig_height = 3.0;
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
    
    hold on;    
    yyaxis right;    
    bar(log10(summary.allele_freqs), 'black', 'BarWidth', 1, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    set(get(gca, 'YAxis'), 'FontSize', 5, 'Color', 'black');
    set(gca, 'ytick', [0:ceil(max(log10(summary.allele_freqs)))]);
    ylabel('log_{10} Transcripts');
    box off;
    
    yyaxis left;
    plot([0:length(summary.alleles)]+0.5, [0; cumsum(summary.allele_freqs)/sum(summary.allele_freqs)], 'Color', 'red', 'LineWidth', 1);
    set(get(gca, 'YAxis'), 'FontSize', 5, 'Color', 'black');    
    ylabel('CDF (over transcripts)', 'FontSize', 6);
    
    set(gca, 'SortMethod', 'depth');
    
    set(get(gca, 'XAxis'), 'FontSize', 5);
    xlabel('Ranked Alleles');
    
    title(title_str, 'FontSize', 6, 'FontWeight', 'normal');
  
    hold off;
    axis tight; 
    box off;
    
    set(gca, 'LineWidth', 1.0);
    set(gca, 'TickDir', 'out');
    
    left_margin = 0.7;
    top_margin = 0.3;    
    bottom_margin = 0.7;        
    right_margin = 0.7;        
    
    set(sp, 'Units', 'centimeters', 'Position', [left_margin, bottom_margin, fig_width-left_margin-right_margin fig_height-top_margin-bottom_margin]);

end

