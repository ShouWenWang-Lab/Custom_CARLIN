


%% construct the object
obj=test_CARLIN_def;

% one way to test whether it works is to call 
%% obj.seq.CARLIN 
% after running this, and compare the sequence with the plasmid 


%% only change the following part


%% % test Tigre 
% A lot of redundant definitions here, but helpful to compute
% them up front once in the singleton definition, rather than
% repeatedly on the fly later.

% Set counts of different element types
obj.N.segments = size(obj.segments,1);
obj.N.pams     = obj.N.segments; % before, obj.N.pams     = obj.N.segments-1; 
obj.N.motifs   = 2*obj.N.segments + obj.N.pams + 2;

% Define types of different motifs
% SW: account for the change of obj.N.pams
obj.motifs.prefix   = 1;
obj.motifs.consites = [2:3:obj.N.motifs-1]; % [2:3:obj.N.motifs]
obj.motifs.cutsites = [3:3:obj.N.motifs-1]; %[3:3:obj.N.motifs]
obj.motifs.pams     = [4:3:obj.N.motifs-1]; 
obj.motifs.postfix  = obj.N.motifs;

%% Define standard sequences
obj.seq.Primer5 = obj.Primer5;
obj.seq.Primer3 = obj.Primer3;
obj.seq.SecondarySequence = obj.SecondarySequence;            

obj.seq.prefix   = obj.prefix;
obj.seq.postfix  = obj.postfix;
obj.seq.pam      = obj.pam;
obj.seq.segments = obj.segments([1:obj.N.segments]);

%SW: account for one more bam-seq; temp(2:2:end-1) may need to be adjusted
temp = cell(obj.N.segments+obj.N.pams+2,1);
temp(1) = cellstr(obj.seq.prefix);
temp(2:2:end-1) = obj.seq.segments; % old temp(2:2:end) = obj.seq.segments;
temp(3:2:end-1) = cellstr(obj.seq.pam); 
temp(end) = cellstr(obj.seq.postfix);
obj.seq.CARLIN = horzcat(temp{:});

%% Precompute widths
obj.width.Primer5 = length(obj.Primer5);
obj.width.Primer3 = length(obj.Primer3);
obj.width.SecondarySequence = length(obj.SecondarySequence);
obj.width.prefix  = length(obj.prefix);
obj.width.segment = length(obj.segments{1});
obj.width.postfix = length(obj.postfix);
obj.width.pam     = length(obj.pam);           
obj.width.CARLIN  = length(obj.seq.CARLIN);
obj.width.cutsite = 7;
obj.width.consite = obj.width.segment-obj.width.cutsite;            
obj.width.min_length = obj.width.prefix+obj.width.consite+obj.width.postfix;

% Define custom sequences
obj.seq.consites = cellfun(@(x) x(1:obj.width.consite), obj.seq.segments, 'un', false);
obj.seq.cutsites = cellfun(@(x) x(obj.width.consite+1:end), obj.seq.segments, 'un', false);

% Define endpoints of different elements
temp = cumsum(cellfun(@length, temp));
temp = [[1; 1+temp(1:end-1)] temp];

%SW: account for one more bam-seq; temp(2:2:end-1) may need to be adjusted
obj.bounds.prefix   = temp(1,:);
obj.bounds.segments = temp(2:2:end-1,:); %old: obj.bounds.segments = temp(2:2:end,:);
obj.bounds.pams     = temp(3:2:end-1,:);
obj.bounds.postfix  = temp(end,:);

%SW: account for one more bam-seq; 
obj.bounds.consites  = [obj.bounds.segments(:,1) obj.bounds.segments(:,2)-obj.width.cutsite];
obj.bounds.cutsites  = [obj.bounds.segments(:,1)+obj.width.segment-obj.width.cutsite obj.bounds.segments(:,2)];            

% Precompute bounds of ordered motifs
obj.bounds.ordered                        = zeros(obj.N.motifs,2);
obj.bounds.ordered(obj.motifs.prefix,:)   = obj.bounds.prefix;
obj.bounds.ordered(obj.motifs.consites,:) = obj.bounds.consites;
obj.bounds.ordered(obj.motifs.cutsites,:) = obj.bounds.cutsites;
obj.bounds.ordered(obj.motifs.pams,:)     = obj.bounds.pams;
obj.bounds.ordered(obj.motifs.postfix,:)  = obj.bounds.postfix;

% Precompute BP membership of all motif types
obj.bps.prefix = obj.bounds.prefix(1):obj.bounds.prefix(2);
obj.bps.consite = arrayfun(@(s,e) s:e, obj.bounds.consites(:,1)', obj.bounds.consites(:,2)', 'un', false);
obj.bps.consite = horzcat(obj.bps.consite{:});
obj.bps.cutsite = arrayfun(@(s,e) s:e, obj.bounds.cutsites(:,1)', obj.bounds.cutsites(:,2)', 'un', false);
obj.bps.cutsite = horzcat(obj.bps.cutsite{:});
obj.bps.pam     = arrayfun(@(s,e) s:e, obj.bounds.pams(:,1)', obj.bounds.pams(:,2)', 'un', false);
obj.bps.pam     = horzcat(obj.bps.pam{:});
obj.bps.postfix = obj.bounds.postfix(1):obj.bounds.postfix(2);

% Define colors for different events used in visualization
obj.color = containers.Map();
obj.color('N') = [255, 255, 255];
obj.color('E') = [255, 220, 220];
obj.color('M') = [0, 255, 0];
obj.color('I') = [0, 0, 255];
obj.color('D') = [255, 0, 0];

% Define transparency for different events used in visualization
obj.alpha.prefix = 0.5;
obj.alpha.consite = 0.7;
obj.alpha.cutsite = 1.0;
obj.alpha.pam = 0.9;
obj.alpha.postfix = 0.5;
obj.alpha.overlay = 0.4;

obj.alpha.CARLIN = zeros(1,obj.width.CARLIN);
obj.alpha.CARLIN(obj.bps.prefix)  = obj.alpha.prefix;
obj.alpha.CARLIN(obj.bps.consite) = obj.alpha.consite;
obj.alpha.CARLIN(obj.bps.cutsite) = obj.alpha.cutsite;
obj.alpha.CARLIN(obj.bps.pam)     = obj.alpha.pam;
obj.alpha.CARLIN(obj.bps.postfix) = obj.alpha.postfix;

% Empirically derived NUC44 alignment score thresholds to determine
% a successful match.
obj.match_score.Primer5   = 17; % has 17bp
obj.match_score.Primer3   = 9; % only has 9bp
obj.match_score.SecondarySequence = 16; % the total seq length is 16

open_penalty = cell(obj.N.motifs,1);
open_penalty(obj.motifs.prefix )  = {10*ones(1, obj.width.prefix)};
temp_SW_postfix=ones(1,obj.width.postfix)+9;
open_penalty(obj.motifs.postfix)  = {temp_SW_postfix};
open_penalty(obj.motifs.cutsites) = {[8:-0.5:6.5 6 6 6]}; % Because we expect the alter- ations to be localized to the cutsites, the minimum value of the penalty was set to occur within 3 bp upstream of the PAM sequences with the penalty increasing linearly moving away from this location. 
open_penalty(obj.motifs.consites) = {[9.4 9.6 9.8 10:-0.125:8.5]}; % I extend it to 16bp now
open_penalty(obj.motifs.pams)     = {[8.6 8.8 9 9.2]};
obj.open_penalty = [10 horzcat(open_penalty{:})];

close_penalty = cell(obj.N.motifs,1);
close_penalty(obj.motifs.prefix )  = {10*ones(1, obj.width.prefix)};
close_penalty(obj.motifs.postfix)  = {temp_SW_postfix};
close_penalty(obj.motifs.cutsites) = {[8:-0.5:6.5 6.5 6.5 6.5]};
close_penalty(obj.motifs.consites) = {[9.4 9.6 9.8 10:-0.125:8.5]};
close_penalty(obj.motifs.pams)     = {[8.6 8.8 9 9.2]};
obj.close_penalty = horzcat(close_penalty{:});

sub_score = nuc44;
obj.sub_score = sub_score(1:4,1:4);
            
            












%% Set counts of different element types
% obj.N.segments = size(obj.segments,1);
% obj.N.pams     = obj.N.segments; % before, obj.N.pams     = obj.N.segments-1; 
% obj.N.motifs   = 2*obj.N.segments + obj.N.pams + 2;
% 
% %% Define types of different motifs
% obj.motifs.prefix   = 1;
% obj.motifs.consites = [4:3:obj.N.motifs]; %before: obj.motifs.consites = [2:3:obj.N.motifs];
% obj.motifs.cutsites = [3:3:obj.N.motifs];% obj.motifs.cutsites = [3:3:obj.N.motifs] 
% obj.motifs.pams     = [2:3:obj.N.motifs-1]; %obj.motifs.pams= [4:3:obj.N.motifs-1]
% obj.motifs.postfix  = obj.N.motifs;
% 
% 
% % Define standard sequences
% obj.seq.Primer5 = obj.Primer5;
% obj.seq.Primer3 = obj.Primer3;
% obj.seq.SecondarySequence = obj.SecondarySequence;            
% 
% obj.seq.prefix   = obj.prefix;
% obj.seq.postfix  = obj.postfix;
% obj.seq.pam      = obj.pam;
% obj.seq.segments = obj.segments([1:obj.N.segments]);
% 
% %% (SW: mirror the structure)
% temp = cell(obj.N.segments+obj.N.pams+2,1);
% temp(1) = cellstr(obj.seq.prefix);
% temp(3:2:end) = obj.seq.segments;% before, temp(2:2:end) = obj.seq.segments;
% temp(2:2:end-1) = cellstr(obj.seq.pam); % temp(3:2:end-1) = cellstr(obj.seq.pam);
% temp(end) = cellstr(obj.seq.postfix);
% obj.seq.CARLIN = horzcat(temp{:});
% 
% 
% % Precompute widths
% obj.width.Primer5 = length(obj.Primer5);
% obj.width.Primer3 = length(obj.Primer3);
% obj.width.SecondarySequence = length(obj.SecondarySequence);
% obj.width.prefix  = length(obj.prefix);
% obj.width.segment = length(obj.segments{1});
% obj.width.postfix = length(obj.postfix);
% obj.width.pam     = length(obj.pam);           
% obj.width.CARLIN  = length(obj.seq.CARLIN);
% obj.width.cutsite = 7;
% obj.width.consite = obj.width.segment-obj.width.cutsite;            
% obj.width.min_length = obj.width.prefix+obj.width.consite+obj.width.postfix;
% 
% %% Define custom sequences (SW: mirror the structure)
% obj.seq.consites = cellfun(@(x) x(obj.width.consite+1:end), obj.seq.segments, 'un', false); % obj.seq.consites = cellfun(@(x) x(1:obj.width.consite), obj.seq.segments, 'un', false);
% obj.seq.cutsites = cellfun(@(x) x(1:obj.width.consite), obj.seq.segments, 'un', false); % obj.seq.cutsites = cellfun(@(x) x(obj.width.consite+1:end), obj.seq.segments, 'un', false);
% 
% % Define endpoints of different elements
% temp = cumsum(cellfun(@length, temp)); % the length for each segment
% temp = [[1; 1+temp(1:end-1)] temp]; %
% 
% obj.bounds.prefix   = temp(1,:);
% obj.bounds.segments = temp(3:2:end,:); % obj.bounds.segments = temp(2:2:end,:);
% obj.bounds.pams     = temp(2:2:end-1,:); % obj.bounds.pams     = temp(3:2:end-1,:);
% obj.bounds.postfix  = temp(end,:);
% 
% %% previous version. Now, switch "cutsite" and "consites"
% %obj.bounds.consites  = [obj.bounds.segments(:,1) obj.bounds.segments(:,2)-obj.width.cutsite];
% %obj.bounds.cutsites  = [obj.bounds.segments(:,1)+obj.width.segment-obj.width.cutsite obj.bounds.segments(:,2)];            
% obj.bounds.cutsites  = [obj.bounds.segments(:,1) obj.bounds.segments(:,2)-obj.width.consite];
% obj.bounds.consites  = [obj.bounds.segments(:,1)+obj.width.segment-obj.width.consite obj.bounds.segments(:,2)];   
% 
% % Precompute bounds of ordered motifs
% obj.bounds.ordered                        = zeros(obj.N.motifs,2);
% obj.bounds.ordered(obj.motifs.prefix,:)   = obj.bounds.prefix;
% obj.bounds.ordered(obj.motifs.consites,:) = obj.bounds.consites;
% obj.bounds.ordered(obj.motifs.cutsites,:) = obj.bounds.cutsites;
% obj.bounds.ordered(obj.motifs.pams,:)     = obj.bounds.pams; 
% obj.bounds.ordered(obj.motifs.postfix,:)  = obj.bounds.postfix;
% 
% % Precompute BP membership of all motif types
% obj.bps.prefix = obj.bounds.prefix(1):obj.bounds.prefix(2);
% obj.bps.consite = arrayfun(@(s,e) s:e, obj.bounds.consites(:,1)', obj.bounds.consites(:,2)', 'un', false);
% obj.bps.consite = horzcat(obj.bps.consite{:});
% obj.bps.cutsite = arrayfun(@(s,e) s:e, obj.bounds.cutsites(:,1)', obj.bounds.cutsites(:,2)', 'un', false);
% obj.bps.cutsite = horzcat(obj.bps.cutsite{:});
% obj.bps.pam     = arrayfun(@(s,e) s:e, obj.bounds.pams(:,1)', obj.bounds.pams(:,2)', 'un', false);
% obj.bps.pam     = horzcat(obj.bps.pam{:});
% obj.bps.postfix = obj.bounds.postfix(1):obj.bounds.postfix(2);
% 
% % Define colors for different events used in visualization
% obj.color = containers.Map();
% obj.color('N') = [255, 255, 255];
% obj.color('E') = [255, 220, 220];
% obj.color('M') = [0, 255, 0];
% obj.color('I') = [0, 0, 255];
% obj.color('D') = [255, 0, 0];
% 
% % Define transparency for different events used in visualization
% obj.alpha.prefix = 0.5;
% obj.alpha.consite = 0.7;
% obj.alpha.cutsite = 1.0;
% obj.alpha.pam = 0.9;
% obj.alpha.postfix = 0.5;
% obj.alpha.overlay = 0.4;
% 
% obj.alpha.CARLIN = zeros(1,obj.width.CARLIN);
% obj.alpha.CARLIN(obj.bps.prefix)  = obj.alpha.prefix;
% obj.alpha.CARLIN(obj.bps.consite) = obj.alpha.consite;
% obj.alpha.CARLIN(obj.bps.cutsite) = obj.alpha.cutsite;
% obj.alpha.CARLIN(obj.bps.pam)     = obj.alpha.pam;
% obj.alpha.CARLIN(obj.bps.postfix) = obj.alpha.postfix;
% 
% % Empirically derived NUC44 alignment score thresholds to determine
% % a successful match.
% obj.match_score.Primer5   = 15;
% obj.match_score.Primer3   = 20;
% obj.match_score.SecondarySequence = 20; % test this please, orignal: 30; we found that 22 is the critical number for our data
% 
% open_penalty = cell(obj.N.motifs,1);
% open_penalty(obj.motifs.prefix )  = {10*ones(1, obj.width.prefix)};
% temp_SW_postfix=ones(1,obj.width.postfix)+9;
% open_penalty(obj.motifs.postfix)  = {temp_SW_postfix}; % edits here, before edition: {[8.6:0.2:10]};
% 
% open_penalty(obj.motifs.cutsites) = {[6 6 6 6.5 7 7.5 8]}; % reversed,  {[8:-0.5:6.5 6 6 6]};
% 
% 
% % I extend it to 17bp from 13bp to account for my shortening the PAM to 3bp
% open_penalty(obj.motifs.consites) =  {[8.5:0.125:10 10 10 10 10]}; %reversed,  {[10:-0.125:8.5]};
% 
% % i shorten it to 3bp from 7bp
% %open_penalty(obj.motifs.pams)     = {[9.8:-0.2:8.6]}; %reversed, {[8.6:0.2:9.8]};
% open_penalty(obj.motifs.pams)     = {[9.8 9.6 9.4]};
% obj.open_penalty = [10 horzcat(open_penalty{:})];
% 
% close_penalty = cell(obj.N.motifs,1);
% close_penalty(obj.motifs.prefix )  = {10*ones(1, obj.width.prefix)};
% close_penalty(obj.motifs.postfix)  = {temp_SW_postfix}; %{[8.6:0.2:10]};
% close_penalty(obj.motifs.cutsites) = {[6.5 6.5 6.5 6.5 7 7.5 8]}; %[8:-0.5:6.5 6.5 6.5 6.5]
% % same here, extend it 
% close_penalty(obj.motifs.consites) = {[8.5:0.125:10 10 10 10 10]}; %{[10:-0.125:8.5]}
% %close_penalty(obj.motifs.pams)     = {[9.8:-0.2:8.6]}; %reversed, {[8.6:0.2:9.8]};
% close_penalty(obj.motifs.pams)     = {[9.8 9.6 9.4]};
% obj.close_penalty = horzcat(close_penalty{:});
% 
% sub_score = nuc44;
% obj.sub_score = sub_score(1:4,1:4);


%% do not change this part
disp('get the CARLIN sequence')
obj.seq.CARLIN 

disp('get the whole sequence, including primers')
% 
%         prefix = 'GC';
%         pam = 'ATGG';
%         postfix = 'A';
%         Primer5 = 'ATGTACAAGTAAAGCGGCC'; %remember to update obj.match_score.Primer5; old: CCTAGCCGGGGATCCTCTAGAGTCGAATGTACAAGTAAAGCGGCC
%         Primer3 = 'TCTAGTTGC';
%         SecondarySequence = 'AGTCTGCTGTGTGCCT';
%         

string(obj.Primer5)+string(obj.seq.CARLIN)+string(obj.SecondarySequence) + string(obj.Primer3)

