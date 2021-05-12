classdef Bank
    
    properties (SetAccess = immutable, GetAccess = public)
        sample_names
        summary
        allele_breakdown_by_sample
        sample_map        
        model
    end
        
    methods (Static)
      
        function prog = getCatchAllPath()            
           [folder, ~, ~] = fileparts(mfilename('fullpath'));
           prog = fileread([folder '/CatchAllPath.txt']);           
           assert(exist(prog, 'file') == 2, 'Specified CatchAll path does not contain the executable');
        end
        
        % Clonal pvalue for allele X (Case 1+2 in manuscript)
        % = P(# of transcripts of allele X > 1 | allele X is seen observed, M observations made)        
        % = P(X > 1 , X > 0, M observations made) / P(X > 0, M observations made)
        % = P(X > 1 , M observations made) / P(X > 0, M observations made)
        % = (1-P(X==0,M)-P(X==1,M))/(1-P(X==0,M))
                
        function p = clonal_pvalue_eqn(rate, M)
            if (nargin == 2)
                rate = rate*M;
            end
            p = (1-exp(-rate)-rate.*exp(-rate))./(1-exp(-rate));
        end
        
        function rate = max_rate_for_sig_level(sig_level, N_obs)            
            syms lambda
            eqn = arrayfun(@(p) Bank.clonal_pvalue_eqn(lambda,N_obs)==p, sig_level, 'un', false);
            rate = cellfun(@(x) double(vpasolve(x, lambda)), eqn, 'un', false);
        end
        
        % Frequency pvalue for allele X (Case 3 in manuscript)
        % = P(# of transcripts of allele X >= c | allele X is observed, M observations made)
        % = (1-sum_{i=0}^{c-1} (P(X==i,M))/(1-P(X==0,M))
        function p = frequency_pvalue_eqn(rate, thresh, M)
            if (nargin == 3)
                rate = rate*M;
            end
            p = poisscdf(thresh-1, rate, 'upper')./poisscdf(0, rate, 'upper');
        end
        
        % Eqn (6-7) of Colwell (2012)
        function [mu, sig] = allele_interpolation_eqn(transcripts_sampled, freqs, freq_counts, alleles_estimated)
            transcripts_obs = sum(freqs.*freq_counts);
            assert(all(transcripts_sampled <= transcripts_obs));         
            mu  = arrayfun(@(a) dot(1-(1-a/transcripts_obs).^freqs, freq_counts), transcripts_sampled);
            sig = sqrt(arrayfun(@(a) dot((1-(1-a/transcripts_obs).^freqs).^2, freq_counts), transcripts_sampled) - mu.^2/alleles_estimated);
        end
        
        % Eqn (12) of Colwell (2012)
        function alleles = allele_extrapolation_eqn(transcripts_sampled, alleles_observed, F0, F1, transcripts_observed)
            assert(all(transcripts_sampled >= transcripts_observed));
            a = transcripts_sampled-transcripts_observed;
            alleles = alleles_observed+F0*(1-exp(-(a./transcripts_observed)*(F1/F0)));
        end
        
        % Eqn (14) of Colwell (2012)
        function transcripts = transcript_extrapolation_eqn(alleles_desired, alleles_observed, alleles_estimated, F0, F1, F2, transcripts_observed)
            assert(all(alleles_desired >= alleles_observed) && all(alleles_desired < alleles_estimated));
            g = alleles_desired/alleles_estimated;
            transcripts = transcripts_observed * ( 1 + F1 / F2 * log(F0./((1-g)*alleles_estimated)) );
        end
        
        function bank = Create(samples, sample_names, outdir)
            bank = Bank(samples, sample_names, outdir);
            save(sprintf('%s/Bank.mat', outdir), 'bank');
        end
        
    end
    
    methods (Access = public)
        
        function obj = Bank(samples, sample_names, outdir)
            
            assert(isa(samples, 'ExperimentSummary'));
            
            assert(size(samples,1) == size(sample_names,1), 'Size of sample_names list should match number of samples');
            
            % 1. Pool samples to make the bank
            
            fprintf('Pooling alleles for bank\n');
            
            obj.sample_names = sample_names;
            [obj.summary, obj.sample_map, obj.allele_breakdown_by_sample] = ExperimentSummary.FromMerge(samples);
            
            if (~exist(outdir, 'dir'))
                mkdir(outdir);               
            end            
            plot_summary(obj.summary, outdir);
            
            % 2. Compute extensive and intensive Poisson rates for observed
            % alleles

            is = strcmp(cellfun(@(x) degap(x.get_seq), obj.summary.alleles, 'un', false), CARLIN_def.getInstance.seq.CARLIN);
            model.transcripts.edited = sum(obj.summary.allele_freqs(~is));
            model.extensive.rates = obj.summary.allele_freqs;
            model.intensive.rates = model.extensive.rates / model.transcripts.edited;
            model.extensive.rates(is) = NaN;
            model.intensive.rates(is) = NaN;
            
            % 3. Tabulate frequency distribution for CatchAll input
            
            fprintf('Fitting abundance distribution with CatchAll...\n');
            
            if (all(is))                
                return;
            end
            
            [model.obs.freqs, ~, model.obs.freq_counts] = find(accumarray(obj.summary.allele_freqs(~is),1));
            
            catchall_subdir = [outdir '/CatchAll'];
            if (~exist(catchall_subdir, 'dir'))
                mkdir(catchall_subdir);                
            end            
            catchall_input_file = sprintf('%s/FrequencyCounts.csv', catchall_subdir);
            
            fid = fopen(catchall_input_file, 'wt');
            for i = 1:length(model.obs.freqs)
                fprintf(fid, '%d,%d\n', model.obs.freqs(i), model.obs.freq_counts(i));
            end
            fclose(fid);
            assert(exist(catchall_input_file, 'file')==2);            
            prog = Bank.getCatchAllPath();            

            if (ispc)
                % Windows version of CatchAll is finnicky about slashes :S
                catchall_input_file = strrep(catchall_input_file, '/', '\');
                catchall_subdir = strrep(catchall_subdir, '/', '\');
            end
            
            % 4. Run CatchAll
            
            system([prog ' ' catchall_input_file ' ' catchall_subdir]);
            
            % 5. Fetch the best fit model name
            
            catchall_output_file = sprintf('%s/FrequencyCounts_BestModelsAnalysis.csv', catchall_subdir);
            
            if (~exist(catchall_output_file, 'file'))
                return;
            end
            
            in = cellfun(@(x) strsplit(x,',', 'CollapseDelimiters', false), ...
                 splitlines(fileread(catchall_output_file)), 'un', false);
            
            best_model_row = find(cellfun(@(x) strcmp(x{1}, 'Best Parm Model'), in));
            while (length(in{best_model_row}) < 3 && ~isempty(in{best_model_row}{1}))
                best_model_row = best_model_row+1;
            end
            if (isempty(in{best_model_row}{1}))
                return;
            end
            model.name = in{best_model_row}{2};            
            model.extensive.cutoff  = str2double(in{best_model_row}{3});
            model.transcripts.cutoff = sum(obj.summary.allele_freqs(obj.summary.allele_freqs <= model.extensive.cutoff));
            
            % 6. Fetch the fitted frequencies from the best fit model
            
            in = cellfun(@(x) strsplit(x,',', 'CollapseDelimiters', false), ...
                 splitlines(fileread(sprintf('%s/FrequencyCounts_BestModelsFits.csv', catchall_subdir))), 'un', false);
            in = in(1:end-1);
            
            model.fit.freqs    = cellfun(@(x) str2double(x{1}), in(2:end));
            Tau = find(model.fit.freqs == model.extensive.cutoff);            
            model.fit.freqs       = model.fit.freqs(1:Tau);            
            model.fit.freq_counts = cellfun(@(x) str2double(x{3}), in(2:Tau+1));
            
            % 7. Grab the parameters and estimate from the best fit model.
            
            in = cellfun(@(x) strsplit(x,',', 'CollapseDelimiters', false), ...
                 splitlines(fileread(sprintf('%s/FrequencyCounts_Analysis.csv', catchall_subdir))), 'un', false);
            
            best_model = in{cellfun(@(x) strcmp(x{1}, model.name) && (str2double(x{2}) == model.extensive.cutoff), in)};
            
            model.alleles.observed   = length(obj.summary.alleles)-1;
            model.alleles.cutoff     = str2double(best_model{3});
            model.alleles.estimated  = floor(str2double(best_model{4}));
            model.alleles.unobserved = model.alleles.estimated-model.alleles.observed;            
            model.alleles.SE         = str2double(best_model{5});
            model.alleles.CI_95_LB   = str2double(best_model{6});
            model.alleles.CI_95_UB   = str2double(best_model{7});
            model.gof.chi2           = str2double(best_model{8});
            model.gof.GOF0           = str2double(best_model{9});
            model.gof.GOF5           = str2double(best_model{10});
            model.gof.AIC            = str2double(best_model{11});
            model.gof.AICc           = str2double(best_model{12});
            
            % The first N components in theta are the exponential
            % distribution parameter for each of the N distributions, and
            % the last N components are the linear combination weighting
            % coefficients for the mixture model summing to 1.
            if (any(strcmp(model.name, {'Poisson'; 'SingleExp'})))
                model.extensive.theta = str2double(best_model{13});
                model.extensive.theta(2) = 1;
            elseif (strcmp(model.name, 'TwoMixedExp'))
                model.extensive.theta = cellfun(@(x) str2double(x), best_model(13:15))';
                model.extensive.theta(4) = 1-model.extensive.theta(3);
            elseif (strcmp(model.name, 'ThreeMixedExp'))
                model.extensive.theta = cellfun(@(x) str2double(x), best_model(13:17))';
                model.extensive.theta(6) = 1-sum(model.extensive.theta(4:5));
            elseif (strcmp(model.name, 'FourMixedExp'))
                model.extensive.theta = cellfun(@(x) str2double(x), best_model(13:19))';
                model.extensive.theta(8) = 1-sum(model.extensive.theta(5:7));
            end
        
            % 8. Rescale to get intensive parameters
            
            % CatchAll gives all parameters as extensive (unitless)
            % quantities. We want to use the per-transcript rates and
            % obtain a different Poisson rate, when testing against a
            % sample with a different number of observations, so do some
            % rescaling here.
            
            N_terms = length(model.extensive.theta)/2;            
            model.intensive.cutoff = model.extensive.cutoff / model.transcripts.edited;
            model.intensive.theta  = model.extensive.theta;
            model.intensive.theta(1:N_terms) = model.intensive.theta(1:N_terms) / model.transcripts.edited;
            
            % 9. Compute useful PDFs and CDFs
            
            if (contains(model.name, 'Exp'))
                
                model.extensive.exp_coeff = model.extensive.theta(1:N_terms);
                model.extensive.lin_coeff = model.extensive.theta(N_terms+1:end);            
                model.intensive.exp_coeff = model.intensive.theta(1:N_terms);
                model.intensive.lin_coeff = model.intensive.theta(N_terms+1:end);
                
                model.extensive.rate_pdf = @(x) sum(model.extensive.lin_coeff .* ...
                                                    exp(-x./model.extensive.exp_coeff) ./ ...
                                                    model.extensive.exp_coeff);
                
                model.extensive.rate_cdf = @(x) sum(model.extensive.lin_coeff .* ...
                                                    (1-exp(-x./model.extensive.exp_coeff)));
             
                model.intensive.rate_pdf = @(x) sum(model.intensive.lin_coeff .* ...
                                                    exp(-x./model.intensive.exp_coeff) ./ ...
                                                    model.intensive.exp_coeff);
                
                model.intensive.rate_cdf = @(x) sum(model.intensive.lin_coeff .* ...
                                                    (1-exp(-x./model.intensive.exp_coeff)));

                model.extensive.count_pdf = @(j) sum(model.extensive.lin_coeff ./ ...
                                                     (1+model.extensive.exp_coeff) .* ...
                                                     (model.extensive.exp_coeff./(1+model.extensive.exp_coeff)).^j);
                                                 
                model.extensive.count_cdf = @(j) sum(model.extensive.lin_coeff .* ...
                                                     (1-(model.extensive.exp_coeff./(1+model.extensive.exp_coeff)).^j));                                                 

                model.intensive.count_pdf = @(j,M) sum(model.intensive.lin_coeff ./ ...
                                                       (1+M*model.intensive.exp_coeff) .* ...
                                                       (M*model.intensive.exp_coeff./(1+M*model.intensive.exp_coeff)).^j);
                                                   
                model.intensive.count_cdf = @(j,M) sum(model.intensive.lin_coeff .* ...
                                                       (1-(M*model.extensive.exp_coeff./(1+M*model.intensive.exp_coeff)).^j));
                                                   
            else
                
                model.extensive.rate_pdf = @(x) dirac(x-model.extensive.theta(1));
                model.extensive.rate_cdf = @(x) heaviside(x-model.extensive.theta(1));
                
                model.intensive.rate_pdf = @(x) dirac(x-model.intensive.theta(1));
                model.intensive.rate_cdf = @(x) heaviside(x-model.intensive.theta(1));
                
                model.extensive.count_pdf = @(k) poisspdf(k, model.extensive.theta(1));
                model.extensive.count_cdf = @(k) poisscdf(k, model.extensive.theta(1));
                
                model.intensive.count_pdf = @(k,M) poisspdf(k, model.intensive.theta(1)*M);
                model.intensive.count_cdf = @(k,M) poisscdf(k, model.intensive.theta(1)*M);
                
            end
            
            % 10. Compute unobserved rate
            
            % P(lambda | X = 0) = P( X = 0 | lambda ) * P(lambda) / P(X = 0)
            %                   = P( X = 0 | lambda ) * P(lambda) / integral [ P( X = 0 | lambda ) * P(lambda) dlambda ]
                                               
            % <lambda | X = 0> = integral [ lambda * P(lambda | X = 0) dlambda ]
            %                  = integral [ lambda * P( X = 0 | lambda ) * P(lambda) dlambda ] /
            %                    integral [          P( X = 0 | lambda ) * P(lambda) dlambda ]
            
            
            model.extensive.mean_unobs_rate = integral(@(x) x.*exp(-x).*model.extensive.rate_pdf(x), 0, max(model.extensive.rates)+1) / ...
                                              integral(@(x)    exp(-x).*model.extensive.rate_pdf(x), 0, max(model.extensive.rates)+1);
            
            model.intensive.mean_unobs_rate = model.extensive.mean_unobs_rate / model.transcripts.edited;
            
            obj.model = model;
        
        end
        
        function [mu, sig] = interpolate_alleles(obj, transcripts_sampled)
            assert(all(transcripts_sampled <= obj.model.transcripts.edited));
            [mu, sig] = Bank.allele_interpolation_eqn(transcripts_sampled, obj.model.obs.freqs, obj.model.obs.freq_counts, obj.model.alleles.estimated);
        end
        
        function [mu, CI_95_LB, CI_95_UB] = extrapolate_alleles(obj, transcripts_sampled)
            assert(all(transcripts_sampled >= obj.model.transcripts.edited));
            F1 =  obj.model.obs.freq_counts(obj.model.obs.freqs==1);
            mu = Bank.allele_extrapolation_eqn(transcripts_sampled, obj.model.alleles.observed, obj.model.alleles.unobserved, ...
                                               F1, obj.model.transcripts.edited);
                                           
            CI_95_LB = Bank.allele_extrapolation_eqn(transcripts_sampled, obj.model.alleles.observed, obj.model.alleles.CI_95_LB-obj.model.alleles.observed, ...
                                                     F1, obj.model.transcripts.edited);
                                                 
            CI_95_UB = Bank.allele_extrapolation_eqn(transcripts_sampled, obj.model.alleles.observed, obj.model.alleles.CI_95_UB-obj.model.alleles.observed, ...
                                                     F1, obj.model.transcripts.edited);            
        end        
       
        function transcripts = extrapolate_transcripts(obj, alleles_desired)
            assert(all(alleles_desired >= obj.model.alleles.observed) && all(alleles_desired < obj.model.alleles.estimated));            
            F1 = obj.model.obs.freq_counts(obj.model.obs.freqs == 1);
            F2 = obj.model.obs.freq_counts(obj.model.obs.freqs == 2);
            transcripts = Bank.transcript_extrapolation_eqn(alleles_desired, obj.model.alleles.observed, obj.model.alleles.estimated, ...
                                                            obj.model.alleles.unobserved, F1, F2, obj.model.transcripts.edited);            
        end
        
        function counts = sample_null_distribution(obj, N_transcripts)
           p = [obj.model.extensive.rates(2:end); obj.model.extensive.mean_unobs_rate*ones(obj.model.alleles.unobserved,1)];
           F = cumsum(p)/sum(p);
           assert(length(F) == obj.model.alleles.estimated);
           counts = accumarray(arrayfun(@(i) find(F >= i, 1, 'first'), rand(N_transcripts,1)), 1, [obj.model.alleles.estimated, 1]);
        end
        
        function [p, rates] = compute_clonal_pvalue(obj, varargin)
            
            if (numel(varargin) == 1)
                summary = varargin{1};
                assert(isa(summary, 'ExperimentSummary'));
                alleles_to_test = summary.alleles;
            elseif (numel(varargin) == 2)
                alleles_to_test = varargin{1};
                Nobs = varargin{2};
            else
                error('Unrecognized arguments when calling allele confidence');
            end            
            
            template_ind = strcmp(cellfun(@(x) degap(x.get_seq), alleles_to_test, 'un', false), CARLIN_def.getInstance.seq.CARLIN);
            
            [is, where] = ismember(cellfun(@(x) degap(x.get_seq), alleles_to_test, 'un', false), ...
                                   cellfun(@(x) degap(x.get_seq), obj.summary.alleles, 'un', false));
            
            rates = zeros(size(is));
            rates(is)  = obj.model.intensive.rates(where(is));
            rates(~is) = obj.model.intensive.mean_unobs_rate;
            
            if (numel(varargin) == 1)            
                Nobs = sum(summary.allele_freqs(~template_ind));            
            end
            p = Bank.clonal_pvalue_eqn(rates, Nobs);
            
            rates(template_ind) = NaN;
            p(template_ind) = 1;
            
        end
        
        function [p, rates] = compute_frequency_pvalue(obj, varargin)
                        
             if (numel(varargin) == 1)
                summary = varargin{1};
                assert(isa(summary, 'ExperimentSummary'));
                alleles_to_test = summary.alleles;
                allele_freqs = summary.allele_freqs;
            elseif (numel(varargin) == 2)
                alleles_to_test = varargin{1};
                allele_freqs = varargin{2};
            else
                error('Unrecognized arguments when calling allele confidence');
            end            
            
            template_ind = strcmp(cellfun(@(x) degap(x.get_seq), alleles_to_test, 'un', false), CARLIN_def.getInstance.seq.CARLIN);
            
            [is, where] = ismember(cellfun(@(x) degap(x.get_seq), alleles_to_test, 'un', false), ...
                                   cellfun(@(x) degap(x.get_seq), obj.summary.alleles, 'un', false));
            
            rates = zeros(size(is));
            rates(is)  = obj.model.intensive.rates(where(is));
            rates(~is) = obj.model.intensive.mean_unobs_rate;
            
            Nobs = sum(allele_freqs(~template_ind));
            
            p = Bank.frequency_pvalue_eqn(rates, allele_freqs, Nobs);
            
            rates(template_ind) = NaN;
            p(template_ind) = 1;
            
        end
        
    end
end