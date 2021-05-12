function [sp, breakdown] = plot_site_decomposition(summary, show_legend, title_str, ylabel_str)

    assert(isa(summary, 'ExperimentSummary'));
    
    left_margin = 0.6;
    right_margin = 0.1;        
    top_margin = 0.3;
    bottom_margin = 0.6;
    
    if (nargin < 4)
        ylabel_str = '';
    end
    
    if (nargin < 3)
        title_str = '';
    end
    
    if (nargin < 2)
        show_legend = false;
    end
    
    if (show_legend)
        plot_height = 2.7;
        fig_width  = 4.275;
    
        leg_height = 1.5;
        leg_width = fig_width;
        fig_height = plot_height+top_margin+bottom_margin+leg_height;
    else
        plot_height = 2.2;
        fig_width  = 3.42;
    
        fig_height = plot_height+top_margin+bottom_margin;
    end
    
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height], ...
           'PaperUnits', 'centimeters', 'PaperSize', [fig_width, fig_height]);        
    sp = subplot(1,1,1);
    
    mut_types = {sprintf('Intersite\nIndel');
                 sprintf('Intersite\nDeletion');
                 sprintf('Intrasite\nIndel');
                 'Insertion'; 
                 sprintf('Intrasite\nDeletion');
                 'Unedited'};
    
    breakdown = zeros(CARLIN_def.getInstance.N.segments, length(mut_types));    
    
    mut_events = cellfun(@(x) Mutation.identify_Cas9_events(x), summary.alleles, 'un', false);
    
    for i = 1:size(summary.alleles,1)
        allele_breakdown = false(size(breakdown));
        if (~isempty(mut_events{i}))
            start_site = arrayfun(@(i) CARLIN_def.coarse_grain_motif(CARLIN_def.find_motif_from_bp(i)), vertcat(mut_events{i}.loc_start));
            end_site   = arrayfun(@(i) CARLIN_def.coarse_grain_motif(CARLIN_def.find_motif_from_bp(i)), vertcat(mut_events{i}.loc_end  ));
            mut_type   = vertcat(mut_events{i}.type);

            % This plot relies on a single event happening at a site, which
            % isn't always true, so pick the most interesting thing.
            % 1) Pick multi-site over single-site
            % 2) Pick compound over insertions over deletion

            priority   = inf(size(mut_type));
            priority(start_site ~= end_site & mut_type == 'C') = 1;
            priority(start_site ~= end_site & mut_type == 'I') = 2;
            priority(start_site ~= end_site & mut_type == 'D') = 3;
            priority(start_site == end_site & (mut_type == 'C' | mut_type == 'M')) = 4;
            priority(start_site == end_site & mut_type == 'I') = 5;
            priority(start_site == end_site & mut_type == 'D') = 6;
            [~, reorder] = sort(priority, 'ascend');
            start_site = start_site(reorder);
            end_site = end_site(reorder);
            mut_type = mut_type(reorder);

            for j = 1:length(mut_type)
                if (any(any(allele_breakdown(start_site(j):end_site(j),:))))
                    continue;
                end
                
                if (start_site(j) ~= end_site(j))
                    if (mut_type(j) == 'C')
                        allele_breakdown(start_site(j):end_site(j), 1) = true;                        
                    elseif (mut_type(j) == 'D')
                        allele_breakdown(start_site(j):end_site(j), 2) = true;
                    end
                else
                    if (mut_type(j) == 'C' || mut_type(j) == 'M')
                        allele_breakdown(start_site(j), 3) = true;
                    elseif (mut_type(j) == 'I')
                        allele_breakdown(start_site(j), 4) = true;
                    elseif (mut_type(j) == 'D')
                        allele_breakdown(start_site(j), 5) = true;
                    end
                end
            end
        end
        allele_breakdown(:,6) = ~any(allele_breakdown, 2);
        
        assert(all(sum(allele_breakdown,2)==1));
        breakdown = breakdown + summary.allele_freqs(i)*allele_breakdown;
    end
    
    assert(all(sum(breakdown,2) == sum(summary.allele_freqs)));
    
    % Purple - indel
    % Red    - deletion
    % Blue   - insertion
    % Alpha1 - intersite
    % Alpha<1 - intrasite
    
    c = [0.5 0 0.5; 1 0 0; 0.5 0 0.5; 0 0 1; 1 0 0; 0.9 0.9 0.9];
    alphaval = [1 1 0.6 1 0.4 1];
    
    % Visually looks nicer if compound, deletion and insertion grouped
    % together
    
    reorder = [1 3 4 2 5 6];
    
    h = bar([1:CARLIN_def.getInstance.N.segments]', breakdown(:,reorder), 'stacked', 'BarWidth', 0.9);
    for i = 1:length(mut_types)
        h(i).FaceColor = c(reorder(i),:);
        h(i).FaceAlpha = alphaval(reorder(i));
    end
    
    axis tight;
    xlim([0.5, CARLIN_def.getInstance.N.segments+0.5]);
    box off;
    
    set(gca, 'xtick', [1:CARLIN_def.getInstance.N.segments], 'xticklabel', num2str([1:CARLIN_def.getInstance.N.segments]'));
    set(get(gca, 'XAxis'), 'FontSize', 5);
    xlabel('CARLIN Target Site', 'FontSize', 6);
    
    ax = gca;
    ax.YAxis.Exponent = floor(log10(sum(summary.allele_freqs)));
    set(get(gca, 'YAxis'), 'FontSize', 5);
    ylabel(ylabel_str, 'FontSize', 6);
        
    title(title_str, 'FontSize', 6, 'FontWeight', 'Normal');
    
    if (show_legend)
        lgd = legend(mut_types(reorder), 'Location', 'SouthOutside', 'NumColumns', 2, 'FontSize', 5);
        legend('boxoff');
    end    
 
    set(sp, 'Units', 'centimeters', 'Position', [left_margin fig_height-top_margin-plot_height fig_width-left_margin-right_margin, plot_height]);        
    
    if (show_legend)
        lgd.Units = 'centimeters';
        lgd.Position = [0 0 leg_width, leg_height];
    end
    
    
end