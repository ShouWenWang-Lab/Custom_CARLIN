function tests = mutation_test
    tests = functiontests(localfunctions);
end

function test_mutation_allele_invariance(testCase)
    
    [folder, ~, ~] = fileparts(mfilename('fullpath'));
    
    golden_mut_list = Mutation.FromFile(sprintf('%s/data/Sanger75Annotations.txt', folder));    
    blacklist = cellfun(@(x) isequal(x, 'X') || isequal(x, '?'), golden_mut_list);    
    alleles = cellfun(@Mutation.apply_mutations, golden_mut_list(~blacklist), 'un', false);
    [~, alleles] = cellfun(@(x) CARLIN_def.cas9_align(degap(x.get_seq)), alleles, 'un', false);
    called_mut_list = cellfun(@Mutation.identify_Cas9_events, alleles, 'un', false);
    verifyEqual(testCase, called_mut_list, golden_mut_list(~blacklist));
        
    Mutation.ToFile(called_mut_list, [folder '/Output'], 'Sanger75Reannotations.txt');    
    golden_strings = splitlines(fileread(sprintf('%s/data/Sanger75Annotations.txt', folder)));
    golden_strings = golden_strings(~blacklist);
    called_strings = splitlines(fileread(sprintf('%s/Output/Sanger75Reannotations.txt', folder)));
    verifyEqual(testCase, called_strings, golden_strings);
   
end