% Modeled after CircularGraph by Paul Kassebaum:
% https://www.mathworks.com/matlabcentral/fileexchange/48576-circulargraph
% Reproducing license below:
%     
% Copyright (c) 2016, The MathWorks, Inc.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the The MathWorks, Inc. nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

classdef (Sealed) plot_stargate < handle
    
    properties (Constant)
        spacer_angle_deg = 2;
        inner_radius = 1;
        outer_radius = 1.1;
        max_width_in_pts = 9;
        min_width_in_pts = 0.35;
        tip_offet_in_deg = 5;
    end
    
    methods (Access = private)
        
        function obj = plot_stargate            
        end
    end

    methods (Static)
        
        % Singleton
        function singleObj = getInstance
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = plot_stargate;
            end
            singleObj = localObj;
        end
        
        function sp = create(summary, Nr, Nc, which_sp)

            assert(isa(summary, 'ExperimentSummary'));

            ref = CARLIN_def.getInstance;
            N_motifs = ref.N.motifs;

            spacer_angle = plot_stargate.spacer_angle_deg;
            const_offset = 90+(360-(ref.width.CARLIN+(N_motifs-1)*spacer_angle))/2;

            val = num2cell(const_offset+ref.bounds.ordered(:,1)+linspace(0,spacer_angle*(N_motifs-1), N_motifs)');
            [plot_params(1:N_motifs).s] = val{:};

            val = num2cell(const_offset+ref.bounds.ordered(:,2)+linspace(0,spacer_angle*(N_motifs-1), N_motifs)');    
            [plot_params(1:N_motifs).e] = val{:};

            plot_params(ref.motifs.prefix).alpha    = ref.alpha.prefix;
            plot_params(ref.motifs.postfix).alpha  = ref.alpha.postfix;

            val = repmat({ref.alpha.consite}, 10, 1);
            [plot_params(ref.motifs.consites).alpha] = val{:};    
            val = repmat({ref.alpha.cutsite}, 10, 1);
            [plot_params(ref.motifs.cutsites).alpha] = val{:};  
            % SW: the update is consistent with the old version
            val = repmat({ref.alpha.pam}, ref.N.pams, 1); % before update: val = repmat({ref.alpha.pam}, 9, 1);
            [plot_params(ref.motifs.pams).alpha]     = val{:};

            indmap = arrayfun(@(i) plot_params(i).s:plot_params(i).e, 1:N_motifs, 'un', false);
            indmap = horzcat(indmap{:});

            inner_radius = plot_stargate.inner_radius;
            outer_radius = plot_stargate.outer_radius;

            if (nargin == 1)
                sp = figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
            else
                sp = subplot(Nr, Nc, which_sp);
            end

            hold on;

            for i = 1:N_motifs        
                theta = linspace(plot_params(i).s, plot_params(i).e, (plot_params(i).e-plot_params(i).s+1)*100);    
                x = [inner_radius * cosd(theta) outer_radius * cosd(fliplr(theta))];
                y = [inner_radius * sind(theta) outer_radius * sind(fliplr(theta))];
                p = patch(x, y, [1,1,1]*plot_params(i).alpha);        
                p.EdgeColor = 'black';
                p.LineWidth = 1;
            end

            tip_x = (inner_radius+outer_radius)/2 * cosd(plot_params(1).s-plot_stargate.tip_offet_in_deg);
            tip_y = (inner_radius+outer_radius)/2 * sind(plot_params(1).s-plot_stargate.tip_offet_in_deg);
            patch([0, 0, tip_x], [inner_radius, outer_radius, tip_y], 'black');

            del_matrix = zeros(ref.width.CARLIN, ref.width.CARLIN);
            ins_matrix = zeros(ref.width.CARLIN, ref.width.CARLIN);
            for i = 1:size(summary.alleles,1)
                mut_events = Mutation.identify_Cas9_events(summary.alleles{i});
                for j = 1:length(mut_events)
                    if (mut_events(j).type=='D' || mut_events(j).type=='C' || mut_events(j).type=='M')
                        del_matrix(mut_events(j).loc_start, mut_events(j).loc_end) = ...
                            del_matrix(mut_events(j).loc_start, mut_events(j).loc_end)+summary.allele_freqs(i);
                    end
                    if (mut_events(j).type=='I' || mut_events(j).type=='C' || mut_events(j).type=='M')
                        L = length(degap(mut_events(j).seq_new));
                        e = min(mut_events(j).loc_start+L-1, ref.width.CARLIN);
                        ins_matrix(mut_events(j).loc_start, e) = ...
                            ins_matrix(mut_events(j).loc_start, e)+summary.allele_freqs(i);
                    end
                end
            end

            [rd, cd, vald] = find(del_matrix);

            [~, idx] = sortrows([cd-rd rd], [1 2], 'descend');
            rd = rd(idx); 
            cd = cd(idx);
            vald = vald(idx);

            [ri, ci, vali] = find(ins_matrix);

            [~, idx] = sortrows([ci-ri ri], [1 2], 'descend');
            ri = ri(idx); 
            ci = ci(idx);
            vali = vali(idx);

            N_del = length(vald);

            r = [rd; ri]; 
            c = [cd; ci];
            val = [vald; vali];
            rel_frac = val/sum(summary.allele_freqs);
            width = log2(rel_frac) + plot_stargate.max_width_in_pts;
            width(width < plot_stargate.min_width_in_pts) = plot_stargate.min_width_in_pts;

            for i = 1:length(val)
                if r(i) ~= c(i)
                    if (abs(indmap(r(i)) - indmap(c(i))) == 180)
                        u = [cosd(indmap(r(i)));cosd(indmap(c(i)))];
                        v = [sind(indmap(r(i)));sind(indmap(c(i)))];                
                    else 
                        u  = [cosd(indmap(r(i)));sind(indmap(r(i)))];
                        v  = [cosd(indmap(c(i)));sind(indmap(c(i)))];
                        x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
                        y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
                        R  = sqrt(x0^2 + y0^2 - 1);
                        thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
                        thetaLim(2) = atan2(v(2)-y0,v(1)-x0);

                        if u(1) >= 0 && v(1) >= 0 
                            theta = [linspace(max(thetaLim),pi,50), linspace(-pi,min(thetaLim),50)].';
                        else
                            theta = linspace(thetaLim(1),thetaLim(2)).';
                        end
                        u = R*cos(theta)+x0;
                        v = R*sin(theta)+y0;
                        u_inc = cosd(indmap(r(i)):indmap(c(i)))';
                        v_inc = sind(indmap(r(i)):indmap(c(i)))';
                        if (~(abs(u_inc(1)-u(end)) < 1e-5 && abs(v_inc(1)-v(end)) < 1e-5))
                            u_inc = flipud(u_inc);
                            v_inc = flipud(v_inc);
                            assert(abs(u_inc(1)-u(end)) < 1e-5 && abs(v_inc(1)-v(end)) < 1e-5);
                        end
                        if (i <= N_del)
                            patch([u; u_inc(2:end-1)], [v; v_inc(2:end-1)], 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.005);
                        end
                    end
                    if (i <= N_del)
                        plot(u, v, '-', 'LineWidth', width(i), 'Color', 'r'); 
                    else
                        plot(u, v, '-', 'LineWidth', width(i), 'Color', 'b'); 
                    end
                else
                    u = cosd(indmap(r(i)));
                    v = sind(indmap(r(i)));
                    if (i <= N_del)
                        scatter(u, v, width(i), 'r', 'filled', 'Marker', 'o'); 
                    else
                        scatter(u, v, width(i), 'b', 'filled', 'Marker', 'o'); 
                    end
                end
            end

            axis equal; 
            axis tight;
            box off;
            axis off;
            hold off;
        end
    end
end
