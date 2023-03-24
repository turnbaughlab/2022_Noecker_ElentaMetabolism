function [curated_model, curated_solution] = standard_elenta_curations(model, met_results)
% Output
%        curated_model: reconstruction with added/edited reactions and bounds
%
%        curated_solution: unconstrained FBA solution with the new
%        reconstruction

%% Metabolite data - get info on metabolites produced or used in any
% condition
if(any(contains(met_results.Properties.VariableNames, 'MaxControl')))
    used_mets_int = unique(met_results((met_results.log2FC < -0.5) & ~strcmp(met_results.MSI, 'NA') & ~strcmp(met_results.MSI, '4') & met_results.MaxControl > 100000,:).BiggInt2);
    prod_mets_int = unique(met_results((met_results.log2FC > 0.5) & ~strcmp(met_results.MSI, 'NA') & ~strcmp(met_results.MSI, '4') & met_results.MaxCulture > 100000,:).BiggInt2);
else
    % Assuming setup of strain dataset, contains Produced and Used
    used_mets_int = unique(met_results(met_results.Used > 0 & ~strcmp(met_results.MSI, 'NA') & ~strcmp(met_results.MSI, '4') ,:).BiggInt2);
    prod_mets_int = unique(met_results(met_results.Produced > 0 & ~strcmp(met_results.MSI, 'NA') & ~strcmp(met_results.MSI, '4') ,:).BiggInt2);
end
used_mets_in_model_int = intersect(used_mets_int, model.mets);
prod_mets_in_model_int = intersect(prod_mets_int, model.mets);
exchange_mets_in_model_int = unique([used_mets_in_model_int; prod_mets_in_model_int]);

% Initial model description
import_comps = model.rxns(printUptakeBound(model));
biomass_rxn = string(checkObjective(model));
BiomassComponentsPos = model.mets(find(model.S(:, strmatch(biomass_rxn, model.rxns)) > 0));
BiomassComponentsNeg = model.mets(find(model.S(:, strmatch(biomass_rxn, model.rxns)) < 0));

%% First, add transporters based on metabolite data
exchange_mets_int_ext = strrep(exchange_mets_in_model_int, '[c]', '[e]');
inBiomassPos = [];
inBiomassNeg = [];
growthBenefit = [];
% Add exchange and transport reactions for compounds depleted from any media
% media - we have evidence E lenta can produce or uptake these
new_sol = optimizeCbModel(model);
for j = 1:length(exchange_mets_int_ext)
    previous_sol = new_sol.f;
    metOrd = findMetIDs(model, exchange_mets_int_ext(j)); 
    metOrdInt = findMetIDs(model, exchange_mets_in_model_int(j));
    inBiomassPos(j) = any(strcmp(BiomassComponentsPos, exchange_mets_in_model_int(j)));
    inBiomassNeg(j) = any(strcmp(BiomassComponentsNeg, exchange_mets_in_model_int(j)));
    if inBiomassPos(j) == 0 % if not in biomass
        if metOrd == 0 %if not in model, add to model
            model = addExchangeRxn(model, exchange_mets_int_ext(j), -1000, 1000);
            met_sub = strrep(exchange_mets_int_ext(j), "[e]", "");
            model = addReaction(model, char(met_sub+'_t2r'), 'reactionFormula', char(string(exchange_mets_in_model_int(j)) + ' -> ' + string(exchange_mets_int_ext(j))), 'lowerBound', -1000, 'upperBound', 1000);
            new_sol = optimizeCbModel(model, 'max');
        else %otherwise, just add transporter
            duplicate = find(sum(model.S ~= 0, 1) == 1 & any(model.S == -1, 1) & any(model.S(metOrd(metOrd~=0), :), 1));
            if length(duplicate)==0
                model = addExchangeRxn(model, exchange_mets_int_ext(j), -1000, 1000);
                new_sol = optimizeCbModel(model, 'max');
                %fva_results = fastFVA(model, 0, 'max', 'ibm_cplex');
            else
                "already can exchange " + exchange_mets_in_model_int(j)
                corresponding_exchange = model.rxns(duplicate);
                corresponding_transport = model.rxns(find(model.S(metOrd, :)));
                corresponding_transport = corresponding_transport(string(corresponding_transport) ~= string(corresponding_exchange));
                % could have used findTrspRxnFromMet(model, metList, compFlag)
                new_sol = optimizeCbModel(model, 'max');
            end
        end
    end
    if new_sol.f < 1
        exchange_mets_int_ext(j)
        break
    end
%    if new_sol.f > 1.01*previous_sol % save info if it increases predicted growth rate
%        growthBenefit = [growthBenefit; exchange_mets_int_ext(j)] ;
%    end
end

%disp("Compounds whose exchange provides a growth benefit: ")
%disp(string(growthBenefit))

%% Additional transport reactions

% Add possible NAD transporter, 3/3/2022
if(~any(ismember(model.rxns, 'EX_nad(e)')))
    model = addExchangeRxn(model, 'nad[e]', -1000, 1000);
    model = addReaction(model, 'nad_t2r', 'reactionFormula', 'nad[c] + h[c] <=> nad[e] + h[e]', 'upperBound', 1000, 'lowerBound', -1000);
    model.rxns{strcmp(model.rxns, 'EX_nad[e]')} = 'EX_nad(e)';
end

% add possible PABA transporter, 3/10/2022
if(~any(ismember(model.rxns, 'EX_4abz(e)')))
    model = addExchangeRxn(model, '4abz[e]', -1000, 1000);
    model = addReaction(model, '4abz_t2r', 'reactionFormula', '4abz[c] + h[c] <=> 4abz[e] + h[e]', 'upperBound', 1000, 'lowerBound', -1000);
    model.rxns{strcmp(model.rxns, 'EX_4abz[e]')} = 'EX_4abz(e)';
end

% Add cysteine transport if not present
if(~any(ismember(model.rxns, 'EX_cys_L(e)')))
    model = addExchangeRxn(model, 'cys_L[e]', -1000, 1000);
    model = addReaction(model, 'CYSt2r', 'reactionFormula', 'cys_L[e] + h[e] <=> cys_L[c] + h[c]', 'subSystem', 'Transport, extracellular', 'upperBound', 1000, 'lowerBound', -1000);
    model.rxns{strcmp(model.rxns, 'EX_cys_L[e]')} = 'EX_cys_L(e)';
end

%% Remove LPS biosynthesis
%lps_rxns = findRxnsFromSubSystem(model, 'Lipopolysaccharide biosynthesis')
lps_rxns = model.rxns(contains(string(model.subSystems), 'Lipopolysaccharide biosynthesis'));
lps_mets = findMetsFromRxns(model, lps_rxns);
bad_mets = intersect(BiomassComponentsNeg, lps_mets);
disp(bad_mets);
bad_mets = bad_mets(~ismember(bad_mets, ["ACP[c]","ctp[c]","h2o[c]","nadp[c]"]));
if(length(bad_mets) > 0)
    model = removeMetabolites(model, string(bad_mets)); %Not removing ACP - just lipidA
end
if(length(lps_rxns) > 0)
    model = removeRxns(model, lps_rxns);
end


%% Arg metabolism refinemnet
if(any(ismember(model.rxns, 'CBMK'))) % if it has CBMKr it doesn't need this
    model = changeRxnBounds(model, 'CBMK', -1000, 'l');
end
% Extra ammonia production allowed
% model = changeRxnBounds(model, 'EX_nh4(e)', 1000, 'u');
% model = changeRxnBounds(model, 'NH4tb', -4000, 'l');
% model = changeRxnBounds(model, 'EX_h(e)', -3000, 'l');
% model = changeRxnBounds(model, 'EX_h2(e)', -3000, 'l');
% model = changeRxnBounds(model, 'H2td', 3000, 'u');

% Reversible NMNAT
%if(any(ismember(model.rxns, 'NMNAT')))
%    model = changeRxnBounds(model, 'NMNAT', -1000, 'l');
%end

% Arginosuccinate and acetate reactions
model.lb(find(strcmp('ARGSSr', model.rxns))) = 0;

% Acetylation of ornithine - metabolomics data
model = addReaction(model, 'ACODAr', 'reactionFormula', 'acorn[c] + h2o[c]  <=> ac[c] + orn[c]', 'subSystem', 'Urea cycle');

% No evidence of ORNTA
model = removeRxns(model, 'ORNTA');
% Remove arginine proton transporter, only use arg-orn exchange
model = removeRxns(model, 'ARGt2');

%% Amino acid biosynthesis (those steps required for growth in EDM)
%Asparagine biosynthesis - g.1738.peg.2298  and 479437.5.peg.1082, 479437.5.2185_12.peg
model = addReaction(model, 'ASNS1', 'reactionFormula', 'asp_L[c] + atp[c] + gln_L[c] + h2o[c]  -> amp[c] + asn_L[c] + glu_L[c] + h[c] + ppi[c] ', 'subSystem', 'Alanine and aspartate metabolism', 'geneRule', 'g.1738.peg.2298');

% Chorismate mutase (phe) - either g.1738.peg.1805 or 479437.5.969_10.peg
model = addReaction(model, 'CHORM', 'reactionFormula', 'chor[c]  -> pphn[c] ', 'subSystem', 'Phenylalanine metabolism', 'geneRule', 'g.1738.peg.1805');

% Alanine biosynthesis
% Elenta 2243 gene is 'g.1738.peg.1224' or "479437.5.1523.peg"
model = addReaction(model, 'ALATA_L', 'reactionFormula', 'akg[c] + ala_L[c] 	<=>	glu_L[c] + pyr[c]', 'subSystem', 'Glutamate metabolism', 'geneRule', 'g.1738.peg.1224');

% Missing lysine biosynthesis steps 
% g.1738.peg.2133, ELEN_RS08135 (Gapmind) -dapA
model = addReaction(model, 'DHDPS', 'reactionFormula', 'aspsa[c] + pyr[c]  -> 23dhdp[c] + 2 h2o[c] + h[c]', 'subSystem', 'Lysine metabolism', 'upperBound', 1000, 'geneRule', 'WP_009307035.1');
% 'g.1738.peg.2001, g.1738.peg.381 - dapB, ELEN_RS02295, ELEN_RS08140 (gm)
model = addReaction(model, 'DHDPRy', 'reactionFormula', '23dhdp[c] + h[c] + nadph[c] 	->	nadp[c] + thdp[c]', 'subSystem', 'Lysine metabolism', 'upperBound', 1000, 'geneRule', 'WP_009306087.1 or WP_009307036.1');
% Make this spontaneous rxn reversible - better than assuming a gap
model.lb(strcmp(model.rxns, 'L2A6ODs')) = -1000;
% 479437.5.788.peg
%model = addReaction(model, 'DAPE', 'reactionFormula', '26dap_LL[c]  <=> 26dap_M[c]', 'subSystem', 'Lysine metabolism', 'lowerBound', -1000, 'upperBound', 1000);
% No gene- this is the gap
%model = addReaction(model, '26DAPLLAT', 'reactionFormula', '26dap_LL[c] + akg[c] 	<=>	glu_L[c] + h2o[c] + h[c] + thdp[c]', 'subSystem', 'Lysine metabolism', 'lowerBound', -1000, 'upperBound', 1000);
%model = remov.lb(contains(model.rxns, 'L2A6ODs')) = 0

% Remove L2A6ODs sink 
model = removeRxns(model, {'DAPabs', 'EX_26dap_M(e)'});

%% Vitamins refinement

% reduce dependency on riboflavin and - 3/10/22 - folate
model.S(find(strcmp(model.mets, 'ribflv[c]')), find(contains(model.rxns, biomass_rxn))) = -0.00005; %-0.0001;
model.S(find(strcmp(model.mets, 'fol[c]')), find(contains(model.rxns, biomass_rxn))) = -0.00001;
model.S(find(strcmp(model.mets, 'cu2[c]')), find(contains(model.rxns, biomass_rxn))) = -0.001;

% also related metabolites - fad and thf
model.S(find(strcmp(model.mets, 'fad[c]')), find(contains(model.rxns, biomass_rxn))) = -0.001; %-0.0001;
model.S(find(strcmp(model.mets, 'thf[c]')), find(contains(model.rxns, biomass_rxn))) = -0.001;

% Glutathione synthesis, this is a gap
model = addReaction(model, 'GLUCYS', 'reactionFormula', 'atp[c] + cys_L[c] + glu_L[c] -> adp[c] + glucys[c] + h[c] + pi[c]', 'subSystem', 'Glutathione metabolism', 'upperBound', 1000, 'lowerBound', -1000);
model = addReaction(model, 'GTHS', 'reactionFormula', 'atp[c] + glucys[c] + gly[c] -> adp[c] + gthrd[c] + h[c] + pi[c]', 'subSystem', 'Glutathione metabolism', 'upperBound', 1000, 'lowerBound', -1000);

% Remove spermidine, doesn't seem to be required or have any relevant
% reactions
model = removeMetabolites(model, 'spmd[c]'); %, 'removeRxnFlag', 'exclusive')
model = removeMetabolites(model, 'spmd[e]'); %, 'removeRxnFlag', 'exclusive')

% Remove zinc, same
model = removeMetabolites(model, 'zn2[c]'); %, 'removeRxnFlag', 'exclusive')

% Remove/constrain CODH reaction, lack of evidence in this condition
model.lb(strcmp(model.rxns, 'CODH_ACS')) = 0;
model.ub(strcmp(model.rxns, 'CODH_ACS')) = 0;


% Add missing fatty acid reactions from 2243 model - this is going to be
% required by all strains
if(~any(ismember(model.rxns, 'EAR180x')))
    model = addReaction(model, 'EAR180x', 'reactionFormula', 'h[c] + nadh[c] + toctd2eACP[c]  -> nad[c] + ocdcaACP[c]', 'subSystem', 'Fatty acid synthesis', 'upperBound', 1000, 'lowerBound', -1000);
    %g.1738.peg.2493
    model = addReaction(model, '3HAD180', 'reactionFormula', '3hoctaACP[c]  -> h2o[c] + toctd2eACP[c]', 'subSystem', 'Fatty acid synthesis', 'upperBound', 1000, 'lowerBound', -1000, 'geneRule', 'g.1738.peg.2493');
    model = addReaction(model, 'EAR160x', 'reactionFormula', 'h[c] + nadh[c] + tpalm2eACP[c]  -> nad[c] + palmACP[c] ', 'subSystem', 'Fatty acid synthesis',  'upperBound', 1000, 'lowerBound', -1000)
    model = addReaction(model, '3HAD160', 'reactionFormula', '3hpalmACP[c]  -> h2o[c] + tpalm2eACP[c]', 'subSystem', 'Fatty acid synthesis', 'upperBound', 1000, 'lowerBound', -1000);
    %g.1738.peg.913
    model = addReaction(model, '3OAS160', 'reactionFormula', 'h[c] + malACP[c] + myrsACP[c] 	->	3opalmACP[c] + ACP[c] + co2[c]', 'subSystem', 'Fatty acid synthesis', 'upperBound', 1000, 'lowerBound', -1000, 'geneRule', 'g.1738.peg.913');
    model = addReaction(model, 'EAR140x', 'reactionFormula', 'h[c] + nadh[c] + tmrs2eACP[c]  -> myrsACP[c] + nad[c] ', 'subSystem', 'Fatty acid synthesis',  'upperBound', 1000, 'lowerBound', -1000)
    model = addReaction(model, '3HAD140', 'reactionFormula', '3hmrsACP[c]  -> h2o[c] + tmrs2eACP[c]', 'subSystem', 'Fatty acid synthesis', 'upperBound', 1000, 'lowerBound', -1000);
    model = addReaction(model, 'EAR120x', 'reactionFormula', 'h[c] + nadh[c] + tddec2eACP[c]  -> ddcaACP[c] + nad[c]', 'subSystem', 'Fatty acid synthesis',  'upperBound', 1000, 'lowerBound', -1000)
    model = addReaction(model, '3HAD120', 'reactionFormula', '3hddecACP[c]  -> h2o[c] + tddec2eACP[c] ', 'subSystem', 'Fatty acid synthesis', 'upperBound', 1000, 'lowerBound', -1000);
    model = addReaction(model, 'EAR100x', 'reactionFormula', 'h[c] + nadh[c] + tdec2eACP[c]  -> dcaACP[c] + nad[c]', 'subSystem', 'Fatty acid synthesis',  'upperBound', 1000, 'lowerBound', -1000);
    model = addReaction(model, '3HAD100', 'reactionFormula', '3hdecACP[c]  -> h2o[c] + tdec2eACP[c]', 'subSystem', 'Fatty acid synthesis', 'upperBound', 1000, 'lowerBound', -1000);
    model = addReaction(model, 'EAR80x', 'reactionFormula', 'h[c] + nadh[c] + toct2eACP[c]  -> nad[c] + ocACP[c]', 'subSystem', 'Fatty acid synthesis',  'upperBound', 1000, 'lowerBound', -1000);
    model = addReaction(model, '3HAD80', 'reactionFormula', '3hoctACP[c]  -> h2o[c] + toct2eACP[c]', 'subSystem', 'Fatty acid synthesis', 'upperBound', 1000, 'lowerBound', -1000);
    model = addReaction(model, 'EAR60x', 'reactionFormula', 'h[c] + nadh[c] + thex2eACP[c]  -> hexACP[c] + nad[c]', 'subSystem', 'Fatty acid synthesis',  'upperBound', 1000, 'lowerBound', -1000);
    model = addReaction(model, '3HAD60', 'reactionFormula', '3hhexACP[c]  -> h2o[c] + thex2eACP[c]', 'subSystem', 'Fatty acid synthesis', 'upperBound', 1000, 'lowerBound', -1000);
    model = addReaction(model, 'EAR40x', 'reactionFormula', 'but2eACP[c] + h[c] + nadh[c]  -> butACP[c] + nad[c]', 'subSystem', 'Fatty acid synthesis',  'upperBound', 1000, 'lowerBound', -1000);
    model = addReaction(model, '3HAD40', 'reactionFormula', '3haACP[c]  -> but2eACP[c] + h2o[c]', 'subSystem', 'Fatty acid synthesis', 'upperBound', 1000, 'lowerBound', -1000);
end

% Change all exchange reactions to use (e) instead of [e]
exchange_rxns = find(contains(model.rxns, "EX_"));
for i=1:length(exchange_rxns)
    model.rxns{exchange_rxns(i)} = char(strrep(model.rxns(exchange_rxns(i)), "[e]", "(e)"));
end

curated_model = model;
curated_solution = optimizeCbModel(curated_model);

end


