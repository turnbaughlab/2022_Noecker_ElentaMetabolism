%% Other type-strain specific curation
function [updated_2243_model] = elenta_2243_curation(model)
%% Apply curations only supported by experimental data and/or genome analysis of type strain

%% Amino acid biosynthesis curation
% Add leucine biosynthesis reaction (most already have it)
if(~any(ismember(model.rxns, 'OMCDC')))
    model = addReaction(model, 'OMCDC', 'reactionFormula', '3c4mop[c] + h[c] -> 4mop[c] + co2[c]', 'subSystem', 'Valine, leucine, and isoleucine metabolism', 'geneRule', 'WP_009607939.1');
end

%Isoleucine
model = addReaction(model, 'CITRAMALS', 'reactionFormula', 'accoa[c] + h2o[c] + pyr[c] -> 2mmal[c] + coa[c] + h[c]', 'subSystem', 'Valine, leucine, and isoleucine metabolism', 'geneRule', '479437.5.399.peg');
% Add Isoleucine biosynthesis (check genes in GapMind)
model = addReaction(model, '2MMALD', 'reactionFormula', '2mmal[c] <=> 2mmale[c] + h2o[c]', 'subSystem', 'Valine, leucine, and isoleucine metabolism', 'geneRule', 'WP_009607939.1 and WP_009306194.1');
model = addReaction(model, '2MMALD2', 'reactionFormula', '2mmale[c] + h2o[c] <=> e3mmal[c]', 'subSystem', 'Valine, leucine, and isoleucine metabolism', 'geneRule', 'WP_009607939.1 and WP_009306194.1');

%% Vitamins

% Missing thiamine biosynthesis reactions in 25559 model: (don't know if these are
% conserved across strains)
if(~any(ismember(model.rxns, 'DXYTST')))
    %thiazole synthase
    % model = addReaction(model, 'DXYTST', 'reactionFormula', 'dxyl5p[c] + imogly[c] + thissh[c] 	<=>	2c4mthzep[c] + 2 h2o[c] + h[c] + this[c]', 'subSystem', 'Thiamine metabolism')
    % thiazole tautomerase - Elen_1745
    model = addReaction(model, 'THZT', 'reactionFormula', '2c4mthzep[c]  <=> cthzp[c]', 'subSystem', 'Thiamine metabolism', 'upperBound', 1000, 'lowerBound', -1000, 'geneRule', 'Elen_1745');
    % alternate thiamine synthesis reaction - Elen_2062
    %model = addReaction(model, 'TMPPP', 'reactionFormula', '2mahmp[c] + cthzp[c] + 2 h[c]  <=> co2[c] + ppi[c] + thmmp[c]', 'subSystem', 'Thiamine metabolism', 'upperBound', 1000, 'lowerBound', -1000, 'geneRule', 'Elen_2062');
    % Use version of TMPPP (Thiamine phosphate diphosphorylase) from the newer VMH database
    model = removeRxns(model, 'TMPPP');
    model = addReaction(model, 'TMPPP', 'reactionFormula', '2mahmp[c] + 4mpetz[c] + h[c]  -> ppi[c] + thmmp[c] ', 'subSystem', 'Thiamine metabolism', 'upperBound', 1000, 'lowerBound', 0, 'geneRule', 'Elen_2062');

    % iminoglycine synthesis  - "479437.5.1802.peg"
    model = addReaction(model, '2IMZS', 'reactionFormula', 'amet[c] + nadph[c] + tyr_L[c]  <=> dad_5[c] + h[c] + imogly[c] + met_L[c] + nadp[c] + pcresol[c] ', 'subSystem', 'Thiamine metabolism', 'upperBound', 1000, 'lowerBound', -1000, 'geneRule', '479437.5.1802.peg');
    % thissh sink
    model = addReaction(model, 'sink_thissh[c]', 'reactionFormula', 'thissh[c] <=> ', 'subSystem', 'Thiamine metabolism', 'upperBound', 1000, 'lowerBound', -1000);
    model = addReaction(model, 'sink_this[c]', 'reactionFormula', 'this[c] <=> ', 'subSystem', 'Thiamine metabolism', 'upperBound', 1000, 'lowerBound', -1000);
    % thiamine phosphatase - Elen_1207 or Elen_3056 from KEGG
    model = addReaction(model, 'THMP', 'reactionFormula', 'h2o[c] + thmmp[c] -> pi[c] + thm[c]', 'subSystem', 'Thiamine metabolism', 'upperBound', 1000, 'lowerBound', -1000, 'geneRule', 'Elen_1207 or Elen_3056');
    % possibility to export excess intermediates
    model = addReaction(model, 'gcald_t2r', 'reactionFormula', 'gcald[e] <=> gcald[c]', 'subSystem', 'Transport, extracellular', 'upperBound', 1000, 'lowerBound', -1000);
    model = addExchangeRxn(model, 'gcald[e]', -1000, 1000);
    model = addReaction(model, '2ahbut_t2r', 'reactionFormula', '2ahbut[c] <=> 2ahbut[e]', 'subSystem', 'Transport, extracellular', 'upperBound', 1000, 'lowerBound', -1000);
    model = addExchangeRxn(model, '2ahbut[e]', -1000, 1000);
    model = addReaction(model, '4hba_t2r', 'reactionFormula', '4hba[c] <=> 4hba[e]', 'subSystem', 'Transport, extracellular', 'upperBound', 1000, 'lowerBound', -1000);
    model = addExchangeRxn(model, '4hba[e]', -1000, 1000);
end

%% Other
% Add putrescine carbamoyltransferase
model = addReaction(model, 'PCBT', 'reactionFormula', 'ncptrc[c] + h[c] + pi[c] <=> ptrc[c] + cbp[c]', 'subSystem', 'Arginine and proline metabolism', 'geneRule', 'WP_009306489.1 or WP_015760997.1');
if(~any(ismember(model.rxns, 'EX_ptrc(e)')))
    model = addExchangeRxn(model, 'ptrc[e]', -1000, 1000);
    model = addReaction(model, 'ptrc_t2r', 'reactionFormula', 'ptrc[c] + h[c] <=> ptrc[e] + h[e]');
    model.rxns{strcmp(model.rxns, 'EX_ptrc[e]')} = 'EX_ptrc(e)';
end

% Better implementation of digoxin rxn - need this to be only for cgr+
% strains
model = addExchangeRxn(model, 'digoxin[e]', -1000, 1000)
model = addExchangeRxn(model, 'dihydro_digoxin[e]', -1000, 1000)
model = addReaction(model, 'DIHYDRO_DIGOXINc', 'reactionFormula', 'digoxin[e] + fadh2[c] -> fad[c] + dihydro_digoxin[e]', 'subSystem', 'Drug metabolism', 'geneRule', 'Elen_2529')

updated_2243_model = model
end
