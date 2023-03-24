%% Clean FBA/FVA results script
initCobraToolbox
solver='ibm_cplex'
%% New cleaned-up version of strain comparison script, using more accurate exchange rate constraints
%met_results = readtable("MediaFiles/biggID_metFCs.txt", 'Format', '%s%s%s%s%s%s%s%s%f%f%f%f%f%s%s')
met_results = readtable('ElentaMediasAll/strainLog2FCs_bigg.txt', 'Format', '%s%s%s%s%s%s%f%f%f%s%s%s')

dir1 = "ElentaFBA/reconstructions/Eggerthella_lenta_2020-5-25/"
egger_models = dir(dir1 + "*.mat")
% Compare with GL experiment results
%media_table = readtable('ElentaMediasAll/medias_mM_FBAw_aaAddIns_2020-9-15.txt')
%media_table = readtable('ElentaMediasAll/medias_mM_FBAfixedDM2_2021-1-21.txt')
media_table = readtable('ElentaMediasAll/medias_mM_FBA2021_04_fattyAcidsVitaminsAgmatine.txt')
% medias_mM_FBA2020-5-19.txt')
% media_table.ExportID(find(strcmp('EX_ser_L(e)', media_table.ExportID))) = {'EX_ser_L[e]'}
n_media = size(media_table, 2)-1
media_names = string(media_table.Properties.VariableNames(2:(n_media+1)))

% Get accurate metabolite data for each strain
    
met_results2243 = met_results(met_results.Name2=="Eggerthella_lenta_DSM2243D" | met_results.Name2=="Eggerthella_lenta_ATCC25559" | met_results.Name2=="Eggerthella_lenta_UCSF2243",:)

%% Read in all models and apply curations
model_list = []
model_rxns = []
model_subsystems = []
for k = 1:length(egger_models)    
    modelfile1 = egger_models(k) %%"Eggerthella_lenta_AGORA_reconstructions/Eggerthella_lenta_RC46F.mat"
    model = readCbModel(char(dir1 + char(modelfile1.name)));
    species_name = strrep(strrep(modelfile1.name,'.mat',''), 'DSM_1', 'DSM1');
    sol1 = optimizeCbModel(model)
    met_results1 = met_results(contains(met_results.Name2, species_name),:)
    if(size(met_results1, 1) == 0)
        if(contains(species_name, '16A'))
            met_results1 = met_results(contains(met_results.Name2, 'Eggerthella_lenta_14A'),:)
        else
            met_results1 = met_results2243
        end
    end
    new_model = standard_elenta_curations(model, met_results1);
    % Additional type-strain-specific curation
    if (k == 14 | k == 22)
        new_model = elenta_2243_curation(new_model);
    end
    model_list = [model_list, new_model];
    new_model = model_list(k);
    model_rxns = [model_rxns; new_model.rxns];
    model_subsystems = [model_subsystems; string(new_model.subSystems)];
    writeCbModel(model_list(k), char(string(species_name)+'.mat'));
end

% Merge 2243 info into 25559
model_list(14) = transfer_reactions(model_list(14), model_list(22));
writeCbModel(model_list(14), 'ElentaMediaMetabolomics/repo/iEL2243_2_items/ElentaATCC25559_curated2021-08.mat')
%writeCbModel(model_list(14), 'ElentaFBA/ElentaATCC25559_curated2021-04.mat')


%% Convert media constraints to appropriate rates
%cellConc = 1.5e7; % we can suppose this is earlier in the growth curve, this is approx an OD of 0.01
cellConc = 1e7; % early-exponential (OD 0.01)
%cellWeight = 0.28*1e-12;
cellWeight = 3.3*1e-13;
t = 40; % Average rate over the whole exponential phase
% See Benchling 3/10/2022 for these values

%current_inf = 1000;
%set_inf = 500;
media_rates = table2array(media_table(:,2:size(media_table, 2)))
for j = 1:length(media_table.ExportID)
    for k = 1:size(media_rates, 2)
        if(media_rates(j,k) ~=0 )
            newRate = conc2Rate(media_rates(j,k), cellConc, t, cellWeight);
            media_rates(j,k) = newRate;
        end
    end
end

% Get list of things present in control media from metabolite data
ctrl_media = readtable('ElentaMediaMetabolomics/presentInCtrls_StrainsGLU.txt')
additional_comps = unique(ctrl_media(~strcmp(ctrl_media.MSI, 'NA') & ctrl_media.MedianCtrlValue > 100000,:).BiggExt2)


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

% Get fluxes for dm1b and dm1bAc0 conditions
j=57
k=14
media_results1 = struct()
media_tab = table(media_table.ExportID) %table(media_table.ExportID(find(media_table.(media_n) > 0)), media_table(find(media_table.(media_n) > 0,media_n)))
media_tab.Var2 = media_rates(:, j)
media_tab.Properties.VariableNames = {'ExportID', 'mM'}
comp_list = media_tab.ExportID(media_tab.mM > 0,:);
additional_comps_mod = intersect(additional_comps, model_list(k).mets);
if(media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) < 0.1 & media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) ~= 0)
  media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) = 0.1;
end

biomass_rxn = checkObjective(model_list(k));
media_tab1 = media_tab

[constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(model_list(k), media_tab1, [], additional_comps_mod);

constrained_fluxes(contains(constrained_fluxes.Var2, 'Lysine metabolism'),:)
constrained_fluxes(contains(constrained_fluxes.Var1, findRxnsFromMets(constrained_mod, '26dap_LL[c]')),:)
constrained_fluxes(contains(constrained_fluxes.Var1, findRxnsFromMets(constrained_mod, 'thdp[c]')),:)
constrained_fluxes(contains(constrained_fluxes.Var1, findRxnsFromMets(constrained_mod, 'lys_L[e]')),:)

% Test range of riboflavin requirements, range of thiamine requirements
ribflv_values = -1e-2:5e-6:-1e-6
growth_rate_results = [];
test_mod = model_list(14);
for val=1:length(ribflv_values)
    val2 = ribflv_values(val)
    test_mod.S(find(strcmp(test_mod.mets, 'ribflv[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    test_mod.S(find(strcmp(test_mod.mets, 'fad[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod, media_tab1, [], additional_comps_mod);
    gr_result = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)));
    growth_rate_results = [growth_rate_results; gr_result];
end
ribflv_table = table(transpose(ribflv_values), growth_rate_results)
writetable(ribflv_table, 'ElentaFBA/ribflv_biomass_vary.txt')


ribflv_values = -1e-2:5e-6:-1e-6
growth_rate_results = [];
test_mod = model_list(14);
for val=1:length(ribflv_values)
    val2 = ribflv_values(val)
    test_mod.S(find(strcmp(test_mod.mets, 'ribflv[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    %test_mod.S(find(strcmp(test_mod.mets, 'fad[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod, media_tab1, [], additional_comps_mod);
    gr_result = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)));
    growth_rate_results = [growth_rate_results; gr_result];
end
ribflv_alone_table = table(transpose(ribflv_values), growth_rate_results)
writetable(ribflv_alone_table, 'ElentaFBA/ribflv_alone_biomass_vary.txt')

ribflv_values = -1e-2:5e-6:-1e-6
growth_rate_results = [];
test_mod = model_list(14);
for val=1:length(ribflv_values)
    val2 = ribflv_values(val)
    %test_mod.S(find(strcmp(test_mod.mets, 'ribflv[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    test_mod.S(find(strcmp(test_mod.mets, 'fad[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod, media_tab1, [], additional_comps_mod);
    gr_result = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)));
    growth_rate_results = [growth_rate_results; gr_result];
end
fad_alone_table = table(transpose(ribflv_values), growth_rate_results)
writetable(fad_alone_table, 'ElentaFBA/fad_biomass_vary.txt')

% Thiamine - 0.0030965 thmpp[c]
thm_values = -2e-2:1e-4:-1e-4
growth_rate_results = [];
test_mod = model_list(14);
for val=1:length(thm_values)
    val2 = thm_values(val)
    test_mod.S(find(strcmp(test_mod.mets, 'thmpp[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod, media_tab1, [], additional_comps_mod);
    gr_result = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)));
    growth_rate_results = [growth_rate_results; gr_result];
end
thm_table = table(transpose(thm_values), growth_rate_results)

% CoA - 0.003
coa_values = -2e-2:1e-4:-1e-4
growth_rate_results = [];
test_mod = model_list(14);
for val=1:length(coa_values)
    val2 = coa_values(val)
    test_mod.S(find(strcmp(test_mod.mets, 'coa[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod, media_tab1, [], additional_comps_mod);
    gr_result = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)));
    growth_rate_results = [growth_rate_results; gr_result];
end
coa_table = table(transpose(coa_values), growth_rate_results)

% adocbl - 0.003
cbl_values = -2e-2:1e-4:-1e-4
growth_rate_results = [];
test_mod = model_list(14);
for val=171:length(cbl_values)
    val2 = cbl_values(val)
    test_mod.S(find(strcmp(test_mod.mets, 'adocbl[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod, media_tab1, [], additional_comps_mod);
    if class(constrained_sol) == 'double'
        gr_result = 0
    else
        gr_result = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)));
    end
    growth_rate_results = [growth_rate_results; gr_result];
end
cbl_table = table(transpose(cbl_values), growth_rate_results)
% Zeros here are because of byproduct of cbl on the opposite side of the
% biomass equation

% Folate related compounds
fol_values = -1e-4:1e-6:-1e-6 %-2e-2:1e-4:-1e-4
growth_rate_results = [];
test_mod = model_list(14);
for val=1:length(fol_values)
    val2 = fol_values(val)
    %test_mod.S(find(strcmp(test_mod.mets, '10fthf[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    %test_mod.S(find(strcmp(test_mod.mets, '5mthf[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    test_mod.S(find(strcmp(test_mod.mets, 'fol[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;    
    [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod, media_tab1, [], additional_comps_mod);
    if class(constrained_sol) == 'double'
        gr_result = 0
    else
        gr_result = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)));
    end
    growth_rate_results = [growth_rate_results; gr_result];
end
fol_table = table(transpose(fol_values), growth_rate_results)

atp_values = -100:4:-40
growth_rate_results = [];
test_mod = model_list(14);
for val=1:length(atp_values)
    val2 = atp_values(val)
    %test_mod.S(find(strcmp(test_mod.mets, '10fthf[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    %test_mod.S(find(strcmp(test_mod.mets, '5mthf[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    test_mod.S(find(strcmp(test_mod.mets, 'atp[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;    
    [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod, media_tab1, [], additional_comps_mod);
    if class(constrained_sol) == 'double'
        gr_result = 0
    else
        gr_result = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)));
    end
    growth_rate_results = [growth_rate_results; gr_result];
end
atp_table = table(transpose(atp_values), growth_rate_results)
writetable(atp_table, 'ElentaFBA/atp_biomass_vary.txt')


biomass_rxn = 'bio1'
all_biomass_components_neg = model_list(14).mets(find(model_list(14).S(:, strmatch(biomass_rxn, model_list(14).rxns)) < 0));
all_biomass_coefs_neg = model_list(14).S(find(model_list(14).S(:, strmatch(biomass_rxn, model_list(14).rxns)) < 0), strmatch(biomass_rxn, model_list(14).rxns));
biomass_double = 2*all_biomass_coefs_neg
biomass_10x = 10*all_biomass_coefs_neg
growth_rate_results_2x = [];
growth_rate_results_10x = [];

for val=1:length(all_biomass_components_neg)
    test_mod = model_list(14);
    val2 = biomass_double(val)
    val10 = biomass_10x(val)
    comp = all_biomass_components_neg(val)
    %test_mod.S(find(strcmp(test_mod.mets, '10fthf[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    %test_mod.S(find(strcmp(test_mod.mets, '5mthf[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    test_mod.S(find(strcmp(test_mod.mets, comp)), find(contains(test_mod.rxns, 'bio1'))) = val2;    
    [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod, media_tab1, [], additional_comps_mod);
    if class(constrained_sol) == 'double'
        gr_result = 0
    else
        gr_result = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)));
    end
    growth_rate_results_2x = [growth_rate_results_2x; gr_result];
    % Same thing with 10x
    test_mod.S(find(strcmp(test_mod.mets, comp)), find(contains(test_mod.rxns, 'bio1'))) = val10;    
    [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod, media_tab1, [], additional_comps_mod);
    if class(constrained_sol) == 'double'
        gr_result = 0
    else
        gr_result = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)));
    end
    growth_rate_results_10x = [growth_rate_results_10x; gr_result];
end
results_table = table(all_biomass_components_neg, biomass_double, growth_rate_results_2x, biomass_10x, growth_rate_results_10x);
writetable(results_table, 'ElentaFBA/all_biomass_components_vary.txt')

test_mod2 = model_list(14)
test_mod2.S(find(strcmp(test_mod2.mets, 'atp[c]')), find(contains(test_mod2.rxns, 'bio1'))) = 2*40.1102;
test_mod2.S(find(strcmp(test_mod2.mets, 'adp[c]')), find(contains(test_mod2.rxns, 'bio1'))) = 2*40;
[constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod2, media_tab1, [], additional_comps_mod);


% Folate related compounds
fol_values = -2e-2:1e-4:-1e-4
growth_rate_results = [];
test_mod = model_list(14);
for val=1:length(fol_values)
    val2 = fol_values(val)
    %test_mod.S(find(strcmp(test_mod.mets, '10fthf[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    %test_mod.S(find(strcmp(test_mod.mets, '5mthf[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;
    test_mod.S(find(strcmp(test_mod.mets, '10fthf[c]')), find(contains(test_mod.rxns, 'bio1'))) = val2;    
    [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(test_mod, media_tab1, [], additional_comps_mod);
    if class(constrained_sol) == 'double'
        gr_result = 0
    else
        gr_result = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)));
    end
    growth_rate_results = [growth_rate_results; gr_result];
end
fol_10fthf_table = table(transpose(fol_values), growth_rate_results)
% Do NAD/FAD also

table(string(constrained_sol.rxns(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'arg_L[c]')))), constrained_sol.solution(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'arg_L[c]'))))

table(string(constrained_sol.rxns(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'atp[c]')) & constrained_sol.solution ~= 0)), constrained_sol.solution(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'atp[c]')) & constrained_sol.solution ~= 0), printRxnFormula(constrained_sol, constrained_sol.rxns(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'atp[c]')) & constrained_sol.solution ~= 0)))
59.179/(59.179+0.010999+0.016868+42.222)

table(string(constrained_sol.rxns(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'nad[c]')) & constrained_sol.solution ~= 0)), constrained_sol.solution(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'nad[c]')) & constrained_sol.solution ~= 0), printRxnFormula(constrained_sol, constrained_sol.rxns(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'nad[c]')) & constrained_sol.solution ~= 0)))
table(string(constrained_sol.rxns(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'nadp[c]')) & constrained_sol.solution ~= 0)), constrained_sol.solution(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'nadp[c]')) & constrained_sol.solution ~= 0), printRxnFormula(constrained_sol, constrained_sol.rxns(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'nadp[c]')) & constrained_sol.solution ~= 0)))

table(string(constrained_sol.rxns(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'akg[c]')) & constrained_sol.solution ~= 0)), constrained_sol.solution(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'akg[c]')) & constrained_sol.solution ~= 0), printRxnFormula(constrained_sol, constrained_sol.rxns(contains(constrained_sol.rxns, findRxnsFromMets(model_list(14), 'akg[c]')) & constrained_sol.solution ~= 0)))

constrained_fluxes(contains(constrained_fluxes.Var1, findRxnsFromMets(constrained_mod, 'akg[c]')),:)

constrained_fluxes(contains(constrained_fluxes.Var1, findRxnsFromMets(constrained_mod, 'nadp[c]')),:)
%constrained_fluxes = table(constrained_fluxes, constrained_mod.ub, constrained_mod.lb)
constrained_fluxes(abs(constrained_fluxes.Var3 - constrained_fluxes.Var4) < 0.00001 & constrained_fluxes.Var3 ~= 0 & (constrained_fluxes.Var3-constrained_fluxes.Var5) < 0.01,:)

constrained_fluxes(contains(constrained_fluxes.Var1, findRxnsFromMets(constrained_mod, 'ribflv[c]')),:)
constrained_fluxes(contains(constrained_fluxes.Var1, findRxnsFromMets(constrained_mod, 'fad[c]')),:)
%constrained_fluxes(abs(constrained_fluxes.Var3-constrained_fluxes.Var4) < 1e-6 & constrained_fluxes.Var4 ~= 0,:)
%constrained_fluxes(contains(constrained_fluxes.Var2, 'Thiamine'),:)
%robustnessAnalysis(constrained_mod, 'EX_thm(e)')
writetable(constrained_fluxes, 'ElentaFBA/DM1b_fva_fluxes_25559.txt')
writetable(table(constrained_sol.rxns, constrained_sol.solution), 'ElentaFBA/DM1b_pfba_fluxes_25559.txt')
writetable(table(constrained_mod.rxns, constrained_mod.rxnNames, string(constrained_mod.subSystems), string(constrained_mod.grRules)), 'ElentaFBA/Elenta25559_rxnInfo.txt')
writetable(table(constrained_mod.genes),'ElentaFBA/Elenta25559_geneInfo.txt')
writeCbModel(constrained_mod, 'ElentaFBA/Elenta25559_curated08-22_DM1b_model.mat')
writeSBML(constrained_mod, 'ElentaFBA/Elenta25559_curated08-22_DM1b_model.sbml')
writeSBML(model_list(14), 'ElentaFBA/Elenta25559_curated08-22_unconstrained_model.sbml')


ac_mod = constrained_mod
ac_lim = ac_mod.lb(contains(ac_mod.rxns, 'EX_ac(e)'))
ac_mod.lb(contains(ac_mod.rxns, 'EX_ac(e)')) = 2*ac_lim
[acetate_fluxes, acetate_growthRates] = robustnessAnalysis(ac_mod, 'EX_ac(e)', 50)

arg_mod = constrained_mod
arg_lim = arg_mod.lb(contains(arg_mod.rxns, 'EX_arg_L(e)'))
arg_mod.lb(contains(arg_mod.rxns, 'EX_arg_L(e)')) = 2*arg_lim
[arg_fluxes, arg_growthRates] = robustnessAnalysis(arg_mod, 'EX_arg_L(e)', 50)
constrained_mod.lb(contains(constrained_mod.rxns, 'EX_ac(e)')) = 2*ac_lim
constrained_mod.lb(contains(constrained_mod.rxns, 'EX_arg_L(e)')) = 1.5*arg_lim
[acetate_fluxes, arg_fluxes, arg_ac_growthRates] = doubleRobustnessAnalysis(constrained_mod, 'EX_ac(e)', 'EX_arg_L(e)', 50)
dlmwrite('ElentaFBA/argAcDoubleRobust.txt', arg_ac_growthRates)
writetable(table(acetate_fluxes, arg_fluxes), 'ElentaFBA/argAcDoubleRobust_acArgRates.txt')


j=75
media_results1 = struct()
media_tab = table(media_table.ExportID) %table(media_table.ExportID(find(media_table.(media_n) > 0)), media_table(find(media_table.(media_n) > 0,media_n)))
media_tab.Var2 = media_rates(:, j)
media_tab.Properties.VariableNames = {'ExportID', 'mM'}
comp_list = media_tab.ExportID(media_tab.mM > 0,:);
additional_comps_mod = intersect(additional_comps, model_list(k).mets);
biomass_rxn = checkObjective(model_list(k));
media_tab1 = media_tab
[constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(model_list(k), media_tab1, [], additional_comps_mod);
writetable(constrained_fluxes, 'ElentaFBA/DM1bAc0_fva_fluxes_25559.txt')
writetable(table(constrained_sol.rxns, constrained_sol.solution), 'ElentaFBA/DM1bAc0_pfba_fluxes_25559.txt')

% All knockouts
all_rxns = model_list(14).rxns
all_genes = model_list(14).genes
[grRatio,grRateKO,grRateWT,hasEffect,delRxn,fluxSolution] = singleRxnDeletion(model_list(14))
base_model_ko_results = table(delRxn, grRatio, grRateKO, hasEffect)
writetable(base_model_ko_results, 'ElentaFBA/base_model25559_curated08-22_singleRxnKOresults.txt')
writetable(table(fluxSolution), 'ElentaFBA/base_model25559_curated08-22_singleRxnKOfluxSolutions.txt')

additional_comps_mod = intersect(additional_comps, model_list(14).mets)
j=57
media_results1 = struct()
media_tab = table(media_table.ExportID) %table(media_table.ExportID(find(media_table.(media_n) > 0)), media_table(find(media_table.(media_n) > 0,media_n)))
media_tab.Var2 = media_rates(:, j)
media_tab.Properties.VariableNames = {'ExportID', 'mM'}
comp_list = media_tab.ExportID(media_tab.mM > 0,:);
additional_comps_mod = intersect(additional_comps, model_list(14).mets);
if(media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) < 0.1 & media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) ~= 0)
  media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) = 0.1;
end

[constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(model_list(14), media_tab, [], additional_comps_mod);

[grRatio,grRateKO,grRateWT,hasEffect,delRxn,fluxSolution] = singleRxnDeletion(constrained_mod)
dm1b_ko_results = table(delRxn, grRatio, grRateKO, hasEffect)
writetable(dm1b_ko_results, 'ElentaFBA/dm1b_model25559_curated08-22_singleRxnKOresults.txt')
writetable(table(fluxSolution), 'ElentaFBA/dm1b_model25559_curated08-22_singleRxnKOfluxSolutions.txt')


% Add back CODH activity for comparison
model_codh_ko = model_list(14)
model_codh_ko.lb(contains(model_list(14).rxns, 'CODH_ACS')) = -1000;
model_codh_ko.ub(contains(model_list(14).rxns, 'CODH_ACS')) = 1000;
j=57
media_results1 = struct()
media_tab = table(media_table.ExportID) %table(media_table.ExportID(find(media_table.(media_n) > 0)), media_table(find(media_table.(media_n) > 0,media_n)))
media_tab.Var2 = media_rates(:, j)
media_tab.Properties.VariableNames = {'ExportID', 'mM'}
comp_list = media_tab.ExportID(media_tab.mM > 0,:);
additional_comps_mod = intersect(additional_comps, model_list(14).mets);
biomass_rxn = checkObjective(model_codh_ko);
if(media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) < 0.1 & media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) ~= 0)
  media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) = 0.1;
end
media_tab1 = media_tab
[constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(model_codh_ko, media_tab1, [], additional_comps_mod);
writetable(constrained_fluxes, 'ElentaFBA/DM1b_withCodh_fva_fluxes_25559.txt')
writetable(table(constrained_sol.rxns, constrained_sol.solution), 'ElentaFBA/DM1b_withCodh_pfba_fluxes_25559.txt')

ac_mod = constrained_mod
ac_lim = ac_mod.lb(contains(ac_mod.rxns, 'EX_ac(e)'))
ac_mod.lb(contains(ac_mod.rxns, 'EX_ac(e)')) = 2*ac_lim
[acetate_fluxes_ko, acetate_growthRates_ko] = robustnessAnalysis(ac_mod, 'EX_ac(e)', 50)
arg_mod = constrained_mod
arg_lim = arg_mod.lb(contains(arg_mod.rxns, 'EX_arg_L(e)'))
arg_mod.lb(contains(arg_mod.rxns, 'EX_arg_L(e)')) = 2*arg_lim
[arg_fluxes_ko, arg_growthRates_ko] = robustnessAnalysis(arg_mod, 'EX_arg_L(e)', 50)
constrained_mod.lb(contains(constrained_mod.rxns, 'EX_ac(e)')) = 2*ac_lim
constrained_mod.lb(contains(constrained_mod.rxns, 'EX_arg_L(e)')) = 1.5*arg_lim
[acetate_fluxes2, arg_fluxes2, arg_ac_growthRates] = doubleRobustnessAnalysis(constrained_mod, 'EX_ac(e)', 'EX_arg_L(e)', 50)
dlmwrite('ElentaFBA/argAcDoubleRobust_withCodh.txt', arg_ac_growthRates)
writetable(table(acetate_fluxes2, arg_fluxes2), 'ElentaFBA/argAcDoubleRobust_acArgRates_withCodh.txt')

acetate_robust_table = table(acetate_fluxes, acetate_growthRates, acetate_fluxes_ko, acetate_growthRates_ko)
arg_robust_table = table(arg_fluxes, arg_growthRates, arg_fluxes_ko, arg_growthRates_ko)
writetable(acetate_robust_table, 'ElentaFBA/DM1b_acetateRobustness_withWithoutCODH.txt')
writetable(arg_robust_table, 'ElentaFBA/DM1b_argRobustness_withWithoutCODH.txt')


%% Test agmatine
j=98
media_results1 = struct()
media_tab = table(media_table.ExportID) %table(media_table.ExportID(find(media_table.(media_n) > 0)), media_table(find(media_table.(media_n) > 0,media_n)))
media_tab.Var2 = media_rates(:, j)
media_tab.Properties.VariableNames = {'ExportID', 'mM'}
comp_list = media_tab.ExportID(media_tab.mM > 0,:);
additional_comps_mod = intersect(additional_comps, model_list(k).mets);
additional_comps_mod = additional_comps_mod(~contains(additional_comps_mod, 'acorn[e]'));
additional_comps_mod = additional_comps_mod(~contains(additional_comps_mod, 'arg_L[e]'));
biomass_rxn = checkObjective(model_list(k));
if(media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) < 0.1 & media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) ~= 0)
  media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) = 0.1;
end
media_tab1 = media_tab
[constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(model_list(k), media_tab1, [], additional_comps_mod);
writetable(constrained_fluxes, 'ElentaFBA/DM1bAgm1Arg0_fva_fluxes_25559.txt')
writetable(table(constrained_sol.rxns, constrained_sol.solution), 'ElentaFBA/DM1bAgm1Arg0_pfba_fluxes_25559.txt')

%% Get fluxes in 70-70 RNA-seq condition
j=113
media_results1 = struct()
media_tab = table(media_table.ExportID) %table(media_table.ExportID(find(media_table.(media_n) > 0)), media_table(find(media_table.(media_n) > 0,media_n)))
media_tab.Var2 = media_rates(:, j)
media_tab.Properties.VariableNames = {'ExportID', 'mM'}
comp_list = media_tab.ExportID(media_tab.mM > 0,:);
additional_comps_mod = intersect(additional_comps, model_list(k).mets);
biomass_rxn = checkObjective(model_list(k));
if(media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) < 0.1 & media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) ~= 0)
  media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) = 0.1;
end
media_tab1 = media_tab
[constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(model_list(k), media_tab1, [], additional_comps_mod);
writetable(constrained_fluxes, 'ElentaFBA/DM1b_70arg-70ac_fva_fluxes_25559.txt')
writetable(table(constrained_sol.rxns, constrained_sol.solution), 'ElentaFBA/DM1b_70arg-70ac_pfba_fluxes_25559.txt')


%% Save updates mets, rxns, genes
updated_met_list = model_list(14).mets;
writetable(table(updated_met_list), 'ElentaFBA/ATCC25559curatedMetList.txt')

updated_rxn_list = model_list(14).rxns;
writetable(table(updated_rxn_list), 'ElentaFBA/ATCC25559curatedRxnList.txt')

updated_gene_list = model_list(14).genes;
writetable(table(updated_gene_list), 'ElentaFBA/ATCC25559curatedGeneList.txt')

updated_rxn_info = table(model_list(14).rxns, printRxnFormula(model_list(14), model_list(14).rxns), model_list(14).grRules);
writetable(updated_rxn_info, 'ElentaFBA/ATCC25559curatedRxnInfo.txt')


% Folate test
model_list(14) = addReaction(model_list(14), 'DHFOR2', 'reactionFormula', 'dhf[c] + nadp[c] <=> fol[c] + nadph[c]', 'subSystem', 'Folate metabolism')
model_list(14) = addReaction(model_list(14), 'DHFR', 'reactionFormula', 'dhf[c] + nadph[c] <=> thf[c] + nadp[c]', 'subSystem', 'Folate metabolism')
media_tab1.mM(find(contains(media_tab1.ExportID, 'EX_fol(e)'))) = 0.01
[constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(model_list(k), media_tab1, [], additional_comps_mod);


%% Save models and updated rxn info for all models
for j = 1:length(model_list)
    updated_rxn_info = table(model_list(j).rxns, model_list(j).subSystems, printRxnFormula(model_list(j), model_list(j).rxns), model_list(j).grRules);
    species_name = strrep(model_list(j).description,'.mat','');
    writetable(updated_rxn_info, char("ElentaFBA/updated_strain_mods/"+species_name+"_curatedRxnInfo.txt"));
    writeCbModel(model_list(j), char("ElentaFBA/updated_strain_mods/"+species_name+"_curatedMod.mat"));
    writeSBML(model_list(j), char("ElentaFBA/updated_strain_mods/"+species_name+"_curatedMod.sbml"));
end
