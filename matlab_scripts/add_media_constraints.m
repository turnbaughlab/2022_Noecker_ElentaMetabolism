function [constrained_model, comp_effects, min_growth_solution, solution_flux_table, retained_comps] = add_media_constraints(model, media, prod_limit_mets, additional_comps_mets)
% Output
%        model with import constraints added based on defined media
%
%        min_growth_solution: FBA solution with these constraints with total flux minimized
%        retained_comps: anything that should have been removed but
%        couldn't be

if(any(strcmp(media.ExportID, 'EX_mndn(e)')))
    mndn_rate = media.mM(find(strcmp(media.ExportID, 'EX_mndn(e)')))
    extra_tab = table({'EX_mqn8(e)'}, mndn_rate/2) % Just assumming an arbitrary fraction is converted for now, may modify later
    extra_tab.Properties.VariableNames = {'ExportID', 'mM'}
    media = [media; extra_tab]
end    

gl_comps = unique(string(media.ExportID))
constrained_model = model
import_comps = model.rxns(printUptakeBound(model))
import_constraints = intersect(gl_comps, import_comps)
disp(import_constraints)
dont_constrain = ["EX_biomass(e)"; "EX_h(e)"]

constrain_qual = ["EX_sheme(e)"; "EX_q8(e)"; "EX_gthrd(e)"] %"EX_ptrc(e)"; 
% Also allow things detected by metabolomics in GLsub to be generally
% present
% Only if they are not otherwise specified as being present in the defined
% media 
% Add compounds identified as being present in controls from metabolite
% data to media specification for GLsub

if(length(additional_comps_mets) > 0)
    add_ids = findMetIDs(model, additional_comps_mets)
    for j = 1:length(add_ids)
        exchangeRxn = model.rxns(find(model.S(add_ids(j), : )))
        exchangeRxn = exchangeRxn(strmatch("EX", exchangeRxn))
        if(~any(strcmp(exchangeRxn, media.ExportID)))
            constrain_qual = [constrain_qual; exchangeRxn]
        end
    end
end

import_constraints = unique([import_constraints; constrain_qual])
limiting_minvit = []; %["EX_fe2(e)"] %; "EX_fe3(e)"; "EX_nac(e)";"EX_ribflv(e)"; "EX_thm(e)"; "EX_pnto_R(e)"; "EX_pheme(e)"] "EX_cobalt2(e)"; "EX_cu2(e)"; "EX_fol(e)"; 
import_constraints = unique([import_constraints; limiting_minvit])
% Set compounds not in media to 0 lb
removed_comps = []
new_val = []
retained_comps = []
for j = 1:length(import_comps) % for every imported compound
    if ~any(strcmp(import_comps(j), import_constraints)) %if its max rate is not specified
        if ~any(strcmp(import_comps(j), dont_constrain)) % if its supposed to be constrained
            import_comps(j)
            constrained_model.lb(find(strcmp(model.rxns, import_comps(j)))) = 0
            new_val1 = optimizeCbModel(constrained_model)
            if(new_val1.f > 1e-5)
                new_val = [new_val; new_val1.f]
                removed_comps = [removed_comps; import_comps(j)]
            else
                retained_comps = [retained_comps; import_comps(j)]
                constrained_model.lb(find(strcmp(model.rxns, import_comps(j)))) = -0.1
            end
%        else
%            constrained_model.lb(find(strcmp(model.rxns, import_comps(j)))) = -1000
%            break
        end
    end
end

comp_removal_effects = table(removed_comps, new_val, repmat("Removed", length(new_val), 1))
comp_removal_effects.Properties.VariableNames{'removed_comps'} = 'Comp';
comp_removal_effects.Properties.VariableNames{'new_val'} = 'value';

amt_constrained_mod = constrained_model
new_val_amt = []
constrained_comps = []
% skipped_comps = []
for j = 1:length(import_constraints)
    import_constraints(j)
    if any(strcmp(import_constraints(j), media.ExportID))
        new_max = media(find(strcmp(import_constraints(j), media.ExportID)),:).(2)
    else
        new_max = 1
    end
    if any(strcmp(import_constraints(j), limiting_minvit))
        new_max = 1
    end
    amt_constrained_mod.lb(find(strcmp(model.rxns, import_constraints(j)))) = -1*new_max
    new_val1 = optimizeCbModel(amt_constrained_mod)
    new_val_amt = [new_val_amt; new_val1.f]
    constrained_comps = [constrained_comps; import_constraints(j)]
end

comp_effects = table(constrained_comps, new_val_amt, repmat("Constrained", length(new_val_amt), 1))
comp_effects.Properties.VariableNames{'constrained_comps'} = 'Comp';
comp_effects.Properties.VariableNames{'new_val_amt'} = 'value';

comp_effects = vertcat(comp_removal_effects, comp_effects)

% Limit production of certain metabolites
new_val_prod = []
for i = 1:length(prod_limit_mets)
    exc_rxn = 'EX_' + strrep(string(prod_limit_mets(i)), "[e]", "(e)")
    if any(strcmp(exc_rxn, amt_constrained_mod.rxns))
        amt_constrained_mod.ub(find(strcmp(amt_constrained_mod.rxns, exc_rxn))) = 5
        new_val1 = optimizeCbModel(amt_constrained_mod)
        new_val_prod = [new_val_prod; new_val1.f]
    else 
        new_val_prod = [new_val_prod; -1]
    end
end
comp_prod_effects = table(prod_limit_mets, new_val_prod, repmat("ProducedLimit", length(new_val_prod), 1))
comp_prod_effects.Properties.VariableNames{'prod_limit_mets'} = 'Comp';
comp_prod_effects.Properties.VariableNames{'new_val_prod'} = 'value';
comp_effects = vertcat(comp_effects, comp_prod_effects)

FBAsolution = optimizeCbModel(amt_constrained_mod, 'max');
if(FBAsolution.f > 1e-10) %Avoid numerical issues
    biomass_rxn = string(checkObjective(amt_constrained_mod))
    pfba_amt = changeRxnBounds(amt_constrained_mod, biomass_rxn, FBAsolution.f, 'l');
    [min_flux, min_sol_irrev] = minimizeModelFlux(pfba_amt, 'min')
    min_sol_irrev.solution = min_flux.v
    [minFlux_amt, maxFlux_amt] = fastFVA(pfba_amt,99);
    fva_table_amt = table(amt_constrained_mod.rxns, string(amt_constrained_mod.subSystems), minFlux_amt, maxFlux_amt, amt_constrained_mod.lb, amt_constrained_mod.ub)
else
    min_sol_irrev = 0
    fva_table_amt = 0
end

constrained_model = amt_constrained_mod
min_growth_solution = min_sol_irrev
solution_flux_table = fva_table_amt
% return [amt_constrained_mod, comp_effects, min_sol_irrev, fva_table_amt]
end


