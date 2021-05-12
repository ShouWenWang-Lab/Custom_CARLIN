function fig = plot_highlighted_alleles(summary, N_alleles_to_highlight, Nr, Nc, which_sp)

    assert(isa(summary, 'ExperimentSummary'));

    fig_width = 10.8;
    fig_height = (N_alleles_to_highlight+2)*fig_width/(CARLIN_def.getInstance.width.CARLIN+2);
    
    
    if (nargin == 2)    
        fig = figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], 'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);
        sp = subplot(1,1,1);
        set(sp, 'Units', 'centimeters', 'Position', [0 0 fig_width, fig_height]);
    else
        fig = subplot(Nr, Nc, which_sp);
    end
        
    RGB = get_sequence_coloring(summary.alleles(1+[1:N_alleles_to_highlight]), 'bp');    
    imshow(padarray(RGB, [1 1], 0)); hold on;
    h = imshow(padarray(repmat([CARLIN_def.getInstance.alpha.CARLIN], [size(RGB,1), 1]), [1 1], 0));  hold off;
    set(h, 'AlphaData', CARLIN_def.getInstance.alpha.overlay);
    axis tight;
    box on;
    daspect([fig_height/(N_alleles_to_highlight+2), fig_width/(CARLIN_def.getInstance.width.CARLIN+2), 1]);    
    
end