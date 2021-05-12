function [ratio_by_eventful_UMI,ratio_by_allele,ave_insert_len,ave_del_len]=ins_del_statistics(summary)        

allele_breakdown_by_sample = summary.allele_freqs;
N_samples = size(allele_breakdown_by_sample,2);
mut_events = cellfun(@Mutation.identify_Cas9_events, summary.alleles, 'un', false);
[num_bp_del, num_bp_ins] = cellfun(@(x) arrayfun(@(i) x(i).num_bps_indel, [1:length(x)]'), mut_events, 'un', false);
L_max = max([vertcat(num_bp_del{:}); vertcat(num_bp_ins{:})]);

num_bp_del = cellfun(@(x) nonzeros(x), num_bp_del, 'un', false);
num_bp_ins = cellfun(@(x) nonzeros(x), num_bp_ins, 'un', false);
N_del = cellfun(@(x) length(x), num_bp_del); % number of deletion events, shape (N_allele,)
N_ins = cellfun(@(x) length(x), num_bp_ins); % number of insertion events, shape (N_allele,)
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


ratio_by_eventful_UMI=sum(ins_freq{:})/sum(del_freq{:});  % an allele might correspond to multiple UMI. 
ratio_by_allele=sum(N_ins)/sum(N_del); % an allele might have multiple insertion or deletion events

x1=1:L_max;
y1=horzcat(ins_freq{:});
ave_insert_len=x1*y1/(sum(y1)); %by Eventful UMI

y2=horzcat(del_freq{:});
ave_del_len=x1*y2/(sum(y2)); %by Eventful UMI

