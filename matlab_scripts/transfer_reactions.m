%%% Transfer reactions with genes over from the 2243 model to
%%% the ATCC25559 model if supported and not already there
function [updated_working_model] = transfer_reactions(working_model, source_model)

m2243 = source_model
m25559 = working_model

addit_rxns = setdiff(m2243.rxns, m25559.rxns);
rxns_w_genes = [];
for r = 1:length(addit_rxns)
    foo = findGenesFromRxns(m2243, addit_rxns(r));
    if length(foo{1}) > 0
        rxns_w_genes = [rxns_w_genes; addit_rxns(r)];
    end
end

% These are added already with slightly different naming
rxns_add = rxns_w_genes(~contains(rxns_w_genes, '3HAD'))
% See curations for folate rxn decisions
rxns_add = rxns_add(~contains(rxns_add, '4HTH'))
% Aconitate hydratase implemented as 2 rxns instead of 1
rxns_add = rxns_add(~contains(rxns_add, 'ACONT'))
% Not enough support for N-acetylorn CT vs Orn/putr (not annotated in KEGG
% or Prokka)
rxns_add = rxns_add(~contains(rxns_add, 'ACOCT'))
% Biotin reaction seems plausible but is blocked (no source of btamp so
% will leave it out for now
rxns_add = rxns_add(~contains(rxns_add, 'BACCLi'))
% Not sufficient evidence
rxns_add = rxns_add(~contains(rxns_add, 'CAP_NR'))
rxns_add = rxns_add(~contains(rxns_add, 'CZP_NR'))
% Different name
rxns_add = rxns_add(~contains(rxns_add, 'CBMK'))
% Blocked rxn, leaving out
rxns_add = rxns_add(~contains(rxns_add, 'CYSTS'))
% See curations script for lysine biosynth curation
rxns_add = rxns_add(~contains(rxns_add, 'DAPE'))
% Redundant
rxns_add = rxns_add(~contains(rxns_add, 'ICDHy'))
% See curations for thiamine metabolism rxns supported by genome checks
rxns_add = rxns_add(~contains(rxns_add, 'KARA'))
% These may be more correct in this version but will leave them out to
% avoid redundancy
rxns_add = rxns_add(~contains(rxns_add, 'DMSOR'))
% Blocked reaction, will leave out for now
rxns_add = rxns_add(~contains(rxns_add, 'DHNAOPT'))
% DOPADH not active in 2243, also not sure the cofactors are correct here
rxns_add = rxns_add(~contains(rxns_add, 'DOPADH'))
% See curations for thiamine metabolism rxns supported by genome checks
rxns_add = rxns_add(~contains(rxns_add, 'DXYTST'))
% SEe curations for fatty acid biosynth decisions
rxns_add = rxns_add(~contains(rxns_add, 'ECOAH2'))
% Blocked
rxns_add = rxns_add(~contains(rxns_add, 'GALASE_LACTLe'))
% Has a different name GLUOX in current mod
rxns_add = rxns_add(~contains(rxns_add, 'GLUSy'))
% See curations for thiamine metabolism info
rxns_add = rxns_add(~contains(rxns_add, 'HPPK'))
% Not going to add nitrate reductase reactions because not experimentally supported
% rn
rxns_add = rxns_add(~contains(rxns_add, 'NO'))
% Not well supported
rxns_add = rxns_add(~contains(rxns_add, 'NZP_NR'))
% Not well supported by other annotations
rxns_add = rxns_add(~contains(rxns_add, 'PC'))
% Not well supported by other annotations
rxns_add = rxns_add(~contains(rxns_add, 'PFK(ppi)'))
rxns_add = rxns_add(~contains(rxns_add, 'PFK_3'))
%Blocked
rxns_add = rxns_add(~contains(rxns_add, 'PHYQS'))
% Already has non-reversible version, will leave as is
rxns_add = rxns_add(~contains(rxns_add, 'PIt6b'))
%blocked
rxns_add = rxns_add(~contains(rxns_add, 'RB5PE'))
%blocked
rxns_add = rxns_add(~contains(rxns_add, 'SADT'))
% See curations for thiamine metabolism rxns supported by genome checks
rxns_add = rxns_add(~contains(rxns_add, 'THII'))
rxns_add = rxns_add(~contains(rxns_add, 'TMPPP_1'))
rxns_add = rxns_add(~contains(rxns_add, 'TMPPP'))
% Blocked - however we do see low levels of indoles produced so may be a
% target for future curation
rxns_add = rxns_add(~contains(rxns_add, 'TRPTA'))
% Will stick with 25559 verison of these reactions for now
rxns_add = rxns_add(~contains(rxns_add, 'UPP3MT'))
%Blocked
rxns_add = rxns_add(~contains(rxns_add, 'URFGTT'))
% Not well supported by other annotations
rxns_add = rxns_add(~contains(rxns_add, 'UTPATPT'))
%Already has non-reversible version, will leave for now b/c don't see
%production of xanthine
rxns_add = rxns_add(~contains(rxns_add, 'XANt2'))
% Adding extracellular version of digoxin rxns instead (see curations)
rxns_add = rxns_add(~contains(rxns_add, 'DIHYDRO_DIG'))

% Going to leave out bile acid reactions for now except for the extracellular dehydrogenases - will benefit from
% additional curation for specifics including cytosolic v extracellular
rxns_add = rxns_add(~contains(rxns_add, 'BIAC'))
rxns_add = rxns_add(~contains(rxns_add, 'BICoAL'))
genes_add = []
for r = 1:length(rxns_add)
    genes_add = [genes_add; m2243.grRules(ismember(m2243.rxns, rxns_add(r)))]
end

m_combined = m25559
for r = 1:length(rxns_add)
    formula = printRxnFormula(m2243, rxns_add(r))
    rxnName1 = m2243.rxnNames(ismember(m2243.rxns, rxns_add(r)))
    subsys = m2243.subSystems(ismember(m2243.rxns, rxns_add(r)))
    rxnID = rxns_add(r)
    geneID = genes_add(r)
    geneID = geneID{1}
    m_combined = addReaction(m_combined, rxnID{1}, 'reactionFormula', formula{1}, 'subSystem', subsys{1}, 'reactionName', rxnName1{1}, 'upperBound', 1000, 'lowerBound', -1000, 'geneRule', geneID);
end    

% Add adenylate kinase amp reactions - gene is
% Eggerthella_lenta_DSM2243REF_02197/479437.5.2224_17.peg
% Add CBLAT - gene is g.1738.peg.1381
%
% Add DHORDnad - genes are 479437.5.1718_11.peg, 479437.5.1719_10.peg
% Add DHORDi - genes are g.1738.peg.2265, g.1738.peg.2475
% Add ICDHy - gene is 479437.5.1513_7.pe
% Add IPDPUPT - gene is g.1738.peg.2367
% Add asp amidohydrolase r0127 - gene is 479437.5.795_16.peg
% Add pant kinase r0578 - genes are 479437.5.1780_9.peg and 479437.5.2941_9.peg
% sure, can leave in SUCDimq - g.1738.peg.42, g.1738.peg.611 are listed
% Add ACCOAC, Eggerthella_lenta_DSM2243REF_02015/g.1738.peg.662

updated_working_model = m_combined;
end
