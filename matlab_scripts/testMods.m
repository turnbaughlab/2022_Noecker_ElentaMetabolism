%%% test
dir1="ElentaMediaMetabolomics/repo/iEL2243_2_items/";

test_mod1 = readCbModel(char(dir1+'ElentaATCC25559_curated2021-08.mat'));

test_mod2 = readCbModel(char(dir1+'Elenta25559_curated08-22_DM1b_model.mat'));

test_mod3 = readSBML(char(dir1+'Elenta25559_curated08-22_DM1b_model.sbml.xml'), 1000)

test_mod4 = readSBML(char(dir1+'Elenta25559_curated08-22_unconstrained_model.sbml.xml'), 1000)