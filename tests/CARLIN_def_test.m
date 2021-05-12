function tests = CARLIN_def_test
    tests = functiontests(localfunctions);
end

function test_cutsite_width(testCase)
    ref = CARLIN_def.getInstance;
    verifyEqual(testCase, ref.width.cutsite, 7);
end

function test_CARLIN_length(testCase)
    ref = CARLIN_def.getInstance;
    verifyEqual(testCase, length(ref.seq.CARLIN), 276);
    verifyEqual(testCase, ref.width.CARLIN, 276);
end

function test_start_pos(testCase)
    ref = CARLIN_def.getInstance;
    verifyEqual(testCase, ref.bounds.prefix(1), 1);
end
   
function test_end_pos(testCase)
    ref = CARLIN_def.getInstance;
    verifyEqual(testCase, ref.bounds.postfix(2), ref.width.CARLIN);
end

function test_continuous_indices_segments(testCase)
    ref = CARLIN_def.getInstance;
    actual = cell(ref.N.segments+ref.N.pams+2,1);
    actual{1}         = ref.bounds.prefix;
    actual(2:2:end)   = num2cell(ref.bounds.segments,2);
    actual(3:2:end-1) = num2cell(ref.bounds.pams,2);
    actual{end}       = ref.bounds.postfix;
    actual = cellfun(@(x) x(1):x(2), actual, 'un', false);
    actual = horzcat(actual{:});
    verifyEqual(testCase, actual, [1:ref.width.CARLIN]);
end

function test_bounds_prefix(testCase)
    ref = CARLIN_def.getInstance;
    verifyEqual(testCase, ref.seq.CARLIN(ref.bounds.prefix(1):ref.bounds.prefix(2)), ref.seq.prefix);
end

function test_bounds_consites(testCase)
    ref = CARLIN_def.getInstance;
    expected = cellfun(@(x) x(1:ref.width.consite), ref.seq.segments, 'un', false);    
    actual_from_bounds = arrayfun(@(s,e) ref.seq.CARLIN(s:e), ref.bounds.consites(:,1), ref.bounds.consites(:,2), 'un', false);
    actual_precomputed = ref.seq.consites;
    verifyEqual(testCase, actual_from_bounds, expected);
    verifyEqual(testCase, actual_precomputed, expected);
end

function test_bounds_cutsites(testCase)
    ref = CARLIN_def.getInstance;
    expected = cellfun(@(x) x(ref.width.segment-ref.width.cutsite+1:end), ref.seq.segments, 'un', false);
    actual_from_bounds = arrayfun(@(s,e) ref.seq.CARLIN(s:e), ref.bounds.cutsites(:,1), ref.bounds.cutsites(:,2), 'un', false);
    actual_precomputed = ref.seq.cutsites;
    verifyEqual(testCase, actual_from_bounds, expected);
    verifyEqual(testCase, actual_precomputed, expected);
end

function test_bounds_segments(testCase)
    ref = CARLIN_def.getInstance;   
    actual_from_bounds = arrayfun(@(s,e) ref.seq.CARLIN(s:e), ref.bounds.segments(:,1), ref.bounds.segments(:,2), 'un', false);
    verifyEqual(testCase, actual_from_bounds, ref.seq.segments);
end

function test_bounds_pams(testCase)
    ref = CARLIN_def.getInstance;
    actual_from_bounds = arrayfun(@(s,e) ref.seq.CARLIN(s:e), ref.bounds.pams(:,1), ref.bounds.pams(:,2), 'un', false);
    verifyEqual(testCase, actual_from_bounds, repelem(cellstr(ref.seq.pam), ref.N.pams, 1));
end

function test_bounds_postfix(testCase)
    ref = CARLIN_def.getInstance;
    verifyEqual(testCase, ref.seq.CARLIN(ref.bounds.postfix(1):ref.bounds.postfix(2)), ref.seq.postfix);
end