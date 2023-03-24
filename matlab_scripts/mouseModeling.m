%% Mouse modeling
initCobraToolbox
solver='ibm_cplex'
%met_results = readtable("MediaFiles/biggID_metFCs.txt", 'Format', '%s%s%s%s%s%s%s%s%f%f%f%f%f%s%s')
met_results = readtable('ElentaMediaMetabolomics/repo/mouse_modeling_files/cecal_media_comp.txt', 'Format', '%s%s%s%s%s%f%f')

curated_mod = readCbModel('ElentaMediaMetabolomics/repo/iEL2243_2_items/iEL2243_2.mat')

% Start with just the current transporters
met_results_sub = met_results(string(met_results.ExtInModel) == '1',:)

% for k = 1:length(met_results_sub.BiggIDExt)
 %   curated_mod = addExchangeRxn(curated_mod, exchange_mets_int_ext(j), -1000, 1000);

met_ids = findMetIDs(curated_mod, met_results_sub.BiggIDExt);
exchangeRxns = strings(length(met_ids), 1)
for k = 1:length(met_ids)
    if(met_ids(k) ~= 0)
        foo = curated_mod.rxns(find(curated_mod.S(met_ids(k), : )))
        exchangeRxn = foo(strmatch("EX", foo))
        exchangeRxns(k) = string(exchangeRxn);
    else
        curated_mod = addMetabolite(curated_mod, met_results_sub.BiggIDExt(k))
        curated_mod = addExchangeRxn(curated_mod, met_results_sub.BiggIDExt(k), -1000, 1000);
        met_sub = strrep(met_results_sub.BiggIDExt(k), "[e]", "");
        curated_mod = addReaction(curated_mod, char(met_sub+'_t2r'), 'reactionFormula', char(string( met_results_sub.BiggIDExt(k)) + ' -> ' + string( met_results_sub.BiggIDInt(k))), 'lowerBound', -1000, 'upperBound', 1000);
        new_met_id = findMetIDs(curated_mod, met_results_sub.BiggIDExt(k))
        foo = curated_mod.rxns(find(curated_mod.S(new_met_id, : )))
        exchangeRxn = foo(strmatch("EX", foo))
        exchangeRxns(k) = string(exchangeRxn);
    end
end

% Change all exchange reactions to use (e) instead of [e]
%exchange_rxns = find(contains(curated_mod.rxns, "EX_"));
%for i=1:length(exchange_rxns)
%    curated_mod.rxns{exchange_rxns(i)} = char(strrep(curated_mod.rxns(exchange_rxns(i)), "[e]", "(e)"));
%end

constraint = transpose(repelem(100, length(met_ids)))
for k = 1:length(met_ids)
    if(met_results_sub.V2(k) > 0.5 & met_results_sub.V1(k) < 7)
        constraint(k) = 1
    end
end

ExportID = exchangeRxns
media_tab = table(ExportID, constraint) %string(met_results_sub.BiggIDExt), curated_mod.rxns(find(curated_mod.S(string(met_results_sub.BiggIDExt), : ))))

[constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(curated_mod, media_tab, [], []);
constrained_fluxes(contains(constrained_fluxes.Var1, 'EX'),:)
constrained_fluxes(contains(constrained_fluxes.Var1, 'EX') & contains(string(constrained_fluxes.Var1), exchangeRxns),:)

constrained_fluxes(contains(constrained_fluxes.Var1, 'EX_') & abs(constrained_fluxes.Var4 - constrained_fluxes.Var5) < 1 & abs(constrained_fluxes.Var4-0) > 1e-8,:)

constrained_fluxes(contains(constrained_fluxes.Var1, findRxnsFromMets(curated_mod, 'atp[c]')),:)

constrained_fluxes(contains(constrained_fluxes.Var1, findRxnsFromMets(curated_mod, 'ACP[c]')),:)
constrained_fluxes_combined = table(constrained_fluxes, constrained_sol.solution)

writetable(constrained_fluxes, 'ElentaMediaMetabolomics/repo/mouse_modeling_files/li_fva_fluxes.txt')
writetable(table(constrained_sol.rxns, constrained_sol.solution), 'ElentaMediaMetabolomics/repo/mouse_modeling_files/li_pfba_fluxes.txt')
writetable(table(constrained_mod.rxns, constrained_mod.rxnNames, string(constrained_mod.subSystems), string(constrained_mod.grRules), string(printRxnFormula(curated_mod, constrained_mod.rxns))), 'ElentaMediaMetabolomics/repo/mouse_modeling_files/li_rxnInfo.txt')
writeCbModel(constrained_mod, 'mouse_cecum_curated_constrained2243.mat')




uncurated_mod = readCbModel('Eggerthella_lenta_ATCC_25559.mat')
[constrained_mod2, comp_effects2, constrained_sol2, constrained_fluxes2, retained_comps2] = add_media_constraints(uncurated_mod, media_tab, [], []);
uncurated_mod2 = readCbModel('Eggerthella_lenta_DSM_2243.mat')
[constrained_mod2, comp_effects2, constrained_sol2, constrained_fluxes2, retained_comps2] = add_media_constraints(uncurated_mod2, media_tab, [], []);

% Incorporate other diff mets
met_results_sub = met_results(string(met_results.ExtInModel) == '1' | abs(met_results.V2) > 0.5,:) 

met_ids = findMetIDs(curated_mod, met_results_sub.BiggIDExt);
exchangeRxns = strings(length(met_ids), 1)
for k = 1:length(met_ids)
    if(met_ids(k) ~= 0)
        foo = curated_mod.rxns(find(curated_mod.S(met_ids(k), : )))
        exchangeRxn = foo(strmatch("EX", foo))
        exchangeRxns(k) = string(exchangeRxn);
    else
        curated_mod = addMetabolite(curated_mod, met_results_sub.BiggIDExt(k))
        curated_mod = addExchangeRxn(curated_mod, met_results_sub.BiggIDExt(k), -1000, 1000);
        met_sub = strrep(met_results_sub.BiggIDExt(k), "[e]", "");
        curated_mod = addReaction(curated_mod, char(met_sub+'_t2r'), 'reactionFormula', char(string( met_results_sub.BiggIDExt(k)) + ' -> ' + string( met_results_sub.BiggIDInt(k))), 'lowerBound', -1000, 'upperBound', 1000);
        new_met_id = findMetIDs(curated_mod, met_results_sub.BiggIDExt(k))
        foo = curated_mod.rxns(find(curated_mod.S(new_met_id, : )))
        exchangeRxn = foo(strmatch("EX", foo))
        exchangeRxns(k) = string(exchangeRxn);
    end
end

% OK constrain by metabolite shifts
curated_mod_li = curated_mod
met_ids = findMetIDs(curated_mod, met_results_sub.BiggIDExt);
for k = 1:length(met_ids)
    if(met_results_sub.V2(k) > 0.5)
        med_row = strmatch(exchangeRxns(k), media_tab.ExportID)
        media_tab.constraint(med_row) = 0.1 %% Edit media table
    elseif(met_results_sub.V2(k) < -0.5)
        curated_mod_li.ub(strmatch(exchangeRxns(k), curated_mod_li.rxns)) = 0;
    end
end

[constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(curated_mod_li, media_tab, [], []);
constrained_fluxes(contains(constrained_fluxes.Var1, 'EX'),:)
constrained_fluxes(contains(constrained_fluxes.Var1, 'EX') & abs(constrained_fluxes.Var4 - constrained_fluxes.Var5) < 1,:)


% Test LOO of all DM1 comps in all strains
j=57
media_results1 = struct()
media_tab = table(media_table.ExportID) %table(media_table.ExportID(find(media_table.(media_n) > 0)), media_table(find(media_table.(media_n) > 0,media_n)))
media_tab.Var2 = media_rates(:, j)
media_tab.Properties.VariableNames = {'ExportID', 'mM'}
comp_list = media_tab.ExportID(media_tab.mM > 0,:)
for k=1:length(model_list)
    additional_comps_mod = intersect(additional_comps, model_list(k).mets);
    biomass_rxn = checkObjective(model_list(k));
    comp_outcomes = []
    for i=1:length(comp_list)
        media_tab1 = media_tab;
        media_tab1.mM(contains(media_tab1.ExportID, comp_list(i)),:) = 0;
        % Pretend there is extra folate
        if(media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) < 0.1 & media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) ~= 0)
            media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) = 0.1;
        end
        %media_tab1 = media_tab1(~contains(media_tab1.ExportID, comp_list(i)),:)
        [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(model_list(k), media_tab1, [], additional_comps_mod);
    
        if(class(constrained_sol)=="struct")
            comp_outcomes = [comp_outcomes; constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)))];
        else
            comp_outcomes = [comp_outcomes; 0];
        end
    end
    comp_table = table(comp_list, comp_outcomes)
    mod_id = strrep(model_list(k).description, ".mat", "")
    writetable(comp_table, 'ElentaFBA/DM1b_LOO_predictions_2022-08_'+mod_id+'.txt')

end
