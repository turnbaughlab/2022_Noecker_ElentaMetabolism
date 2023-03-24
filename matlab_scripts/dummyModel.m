%% Run paper analyses with uncurated model

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

met_results2243 = met_results(met_results.Name2=="Eggerthella_lenta_DSM2243D" | met_results.Name2=="Eggerthella_lenta_ATCC25559" | met_results.Name2=="Eggerthella_lenta_UCSF2243",:)

model_list = []
model_rxns = []
model_subsystems = []
% Get all biomass components
% Biomass differences
% BiomassPos = {}
% BiomassNeg = {}
BiomassTables = {}
used_2243 = []
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
            used_2243 = [used_2243, species_name]
        end
    end
    %new_model = standard_elenta_curations(model, met_results1);
    new_model = model
    model_list = [model_list, new_model]
    new_model = model_list(k)
    model_rxns = [model_rxns; new_model.rxns]
    model_subsystems = [model_subsystems; string(new_model.subSystems)]
    biomass_rxn = checkObjective(model_list(k))
    BiomassComponentsPos = model_list(k).mets(find(model_list(k).S(:, strmatch(biomass_rxn, model_list(k).rxns)) > 0))
    BiomassComponentsNeg = model_list(k).mets(find(model_list(k).S(:, strmatch(biomass_rxn, model_list(k).rxns)) < 0))
    BiomassPos{k} = sort(BiomassComponentsPos)
    BiomassNeg{k} = sort(BiomassComponentsNeg)
    BiomassTables{k} = table(model_list(k).mets(find(model_list(k).S(:, strmatch(biomass_rxn, model_list(k).rxns)) ~=0 )), model_list(k).S(find(model_list(k).S(:, strmatch(biomass_rxn, model_list(k).rxns))), strmatch(biomass_rxn, model_list(k).rxns)))
end

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


all_prediction_results = []
k = 14
%for k = 1:length(model_list)
    all_results_mod = []
    media_results = []
    biomass_rxn = checkObjective(model_list(k));
    additional_comps_mod = intersect(additional_comps, model_list(k).mets);

    for j = 1:n_media
        media_results1 = struct();
        media_tab = table(media_table.ExportID); %table(media_table.ExportID(find(media_table.(media_n) > 0)), media_table(find(media_table.(media_n) > 0,media_n)))
        media_tab.Var2 = media_rates(:, j);
        media_tab.Properties.VariableNames = {'ExportID', 'mM'}
        % Pretend there is extra folate
        if(media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) < 0.1 & media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) ~= 0)
            media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) = 0.1;
        end
        [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(model_list(k), media_tab, [], additional_comps_mod);
        if(class(constrained_sol)=="struct")
            media_results(j) = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)))
        else
            media_results(j) = 0
        end
        media_results1.model = constrained_mod
        media_results1.compTable = comp_effects
        media_results1.sol = constrained_sol
        media_results1.fluxTable = constrained_fluxes
        media_results1.retainedComps = retained_comps
        all_results_mod = [all_results_mod, media_results1]
    end
%    all_prediction_results = [all_prediction_results; all_results_mod]
    media_results_table = table(transpose(media_names), transpose(media_results))
    writetable(media_results_table, 'ElentaFBA/experimentalMediasDummyModel_predictions_2022-03_'+strrep(model_list(k).description, ".mat", "")+'.txt')

%end
j=57
k=14
media_results1 = struct()
media_tab = table(media_table.ExportID) %table(media_table.ExportID(find(media_table.(media_n) > 0)), media_table(find(media_table.(media_n) > 0,media_n)))
media_tab.Var2 = media_rates(:, j)
media_tab.Properties.VariableNames = {'ExportID', 'mM'}
comp_list = media_tab.ExportID(media_tab.mM > 0,:)
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
    writetable(comp_table, 'ElentaFBA/DM1b_LOO_predictionsDummyModel_2022-03_'+mod_id+'.txt')

% Repeat with 2243 model

all_prediction_results = []
k = 22
%for k = 1:length(model_list)
    all_results_mod = []
    media_results = []
    biomass_rxn = checkObjective(model_list(k));
    additional_comps_mod = intersect(additional_comps, model_list(k).mets);

    for j = 1:n_media
        media_results1 = struct();
        media_tab = table(media_table.ExportID); %table(media_table.ExportID(find(media_table.(media_n) > 0)), media_table(find(media_table.(media_n) > 0,media_n)))
        media_tab.Var2 = media_rates(:, j);
        media_tab.Properties.VariableNames = {'ExportID', 'mM'}
        % Pretend there is extra folate
        if(media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) < 0.1 & media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) ~= 0)
            media_tab.mM(contains(media_tab.ExportID, 'EX_fol(e)')) = 0.1;
        end
        [constrained_mod, comp_effects, constrained_sol, constrained_fluxes, retained_comps] = add_media_constraints(model_list(k), media_tab, [], additional_comps_mod);
        if(class(constrained_sol)=="struct")
            media_results(j) = constrained_sol.lb(find(strcmp(biomass_rxn, constrained_sol.rxns)))
        else
            media_results(j) = 0
        end
        media_results1.model = constrained_mod
        media_results1.compTable = comp_effects
        media_results1.sol = constrained_sol
        media_results1.fluxTable = constrained_fluxes
        media_results1.retainedComps = retained_comps
        all_results_mod = [all_results_mod, media_results1]
    end
%    all_prediction_results = [all_prediction_results; all_results_mod]
    media_results_table = table(transpose(media_names), transpose(media_results))
    writetable(media_results_table, 'ElentaFBA/experimentalMediasDummyModel_predictions_2022-03_'+strrep(model_list(k).description, ".mat", "")+'.txt')

    all_results_mod(57).compTable %% Thinks it needs citrate
%end

j=57
k=22
media_results1 = struct()
media_tab = table(media_table.ExportID) %table(media_table.ExportID(find(media_table.(media_n) > 0)), media_table(find(media_table.(media_n) > 0,media_n)))
media_tab.Var2 = media_rates(:, j)
media_tab.Properties.VariableNames = {'ExportID', 'mM'}
comp_list = media_tab.ExportID(media_tab.mM > 0,:)
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
    writetable(comp_table, 'ElentaFBA/DM1b_LOO_predictionsDummyModel_2022-03_'+mod_id+'.txt')


% What about carveme or model seed models?
