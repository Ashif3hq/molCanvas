// Literature ¹H / ¹³C shifts for the NMR benchmark (scripts/check-nmr-benchmark.ts).
//
// Each entry lists one row per chemically distinct environment, as [δ ppm, number of nuclei]. When the
// source reports two environments as one line (overlapping multiplets), the line is repeated once per
// environment with its share of the nuclei, so the number of rows is the number of signals a perfect
// equivalence analysis should give. Exchangeable O–H / N–H protons are left out. Diastereotopic protons
// are separate rows. Values are copied as printed in the cited source; nothing is estimated.
//
// Sources (cite by key):
//   F2010  G. R. Fulmer, A. J. M. Miller, N. H. Sherden, H. E. Gottlieb, A. Nudelman, B. M. Stoltz,
//          J. E. Bercaw, K. I. Goldberg, "NMR Chemical Shifts of Trace Impurities: Common Laboratory
//          Solvents, Organics, and Gases in Deuterated Solvents Relevant to the Organometallic Chemist",
//          Organometallics 2010, 29, 2176–2179, Tables 1 and 2 (CDCl₃ column). doi:10.1021/om100106e
//   S2022  T. Stadelmann, C. Balmer, S. Riniker, M.-O. Ebert, "Impact of solvent interactions on ¹H and ¹³C
//          chemical shifts investigated using DFT and a reference dataset recorded in CDCl₃ and CCl₄",
//          Phys. Chem. Chem. Phys. 2022, 24, 23551–23560, Table 1 (CDCl₃; fully assigned, with the two H
//          of every diastereotopic CH₂ listed separately). doi:10.1039/d2cp03205h. Compound numbers below
//          are the paper's; integrations follow from the structures.

export interface BenchEntry {
  name: string;
  smiles: string;
  solvent: 'CDCl3' | 'DMSO-d6' | 'D2O' | 'CD3OD';
  source: string;
  /** What makes it a useful test of equivalence (shown in the report). */
  tests?: string;
  h: [number, number][];
  c: [number, number][];
}

export const SOURCES: Record<string, string> = {
  F2010: 'Fulmer et al., Organometallics 2010, 29, 2176 (doi:10.1021/om100106e), Tables 1–2, CDCl₃',
  S2022: 'Stadelmann et al., Phys. Chem. Chem. Phys. 2022, 24, 23551 (doi:10.1039/d2cp03205h), Table 1, CDCl₃',
};

export const BENCHMARK: BenchEntry[] = [
  // ---- Fulmer 2010 (CDCl₃) ----------------------------------------------------------------------
  { name: 'acetone', smiles: 'CC(C)=O', solvent: 'CDCl3', source: 'F2010', tests: 'homotopic methyls', h: [[2.17, 6]], c: [[207.07, 1], [30.92, 2]] },
  { name: 'acetonitrile', smiles: 'CC#N', solvent: 'CDCl3', source: 'F2010', h: [[2.10, 3]], c: [[116.43, 1], [1.89, 1]] },
  { name: 'acetic acid', smiles: 'CC(=O)O', solvent: 'CDCl3', source: 'F2010', h: [[2.10, 3]], c: [[175.99, 1], [20.81, 1]] },
  { name: 'benzene', smiles: 'c1ccccc1', solvent: 'CDCl3', source: 'F2010', tests: 'six-fold ring symmetry', h: [[7.36, 6]], c: [[128.37, 6]] },
  { name: 'tert-butyl alcohol', smiles: 'CC(C)(C)O', solvent: 'CDCl3', source: 'F2010', h: [[1.28, 9]], c: [[69.15, 1], [31.25, 3]] },
  { name: 'chloroform', smiles: 'ClC(Cl)Cl', solvent: 'CDCl3', source: 'F2010', h: [[7.26, 1]], c: [[77.36, 1]] },
  { name: 'cyclohexane', smiles: 'C1CCCCC1', solvent: 'CDCl3', source: 'F2010', tests: 'ring flip averages axial/equatorial', h: [[1.43, 12]], c: [[26.94, 6]] },
  { name: '1,2-dichloroethane', smiles: 'ClCCCl', solvent: 'CDCl3', source: 'F2010', h: [[3.73, 4]], c: [[43.50, 2]] },
  { name: 'dichloromethane', smiles: 'ClCCl', solvent: 'CDCl3', source: 'F2010', tests: 'homotopic CH₂', h: [[5.30, 2]], c: [[53.52, 1]] },
  { name: 'diethyl ether', smiles: 'CCOCC', solvent: 'CDCl3', source: 'F2010', tests: 'enantiotopic CH₂', h: [[3.48, 4], [1.21, 6]], c: [[65.91, 2], [15.20, 2]] },
  { name: 'diglyme', smiles: 'COCCOCCOC', solvent: 'CDCl3', source: 'F2010', h: [[3.65, 4], [3.57, 4], [3.39, 6]], c: [[71.90, 2], [70.51, 2], [59.01, 2]] },
  { name: 'DMF', smiles: 'CN(C)C=O', solvent: 'CDCl3', source: 'F2010', tests: 'slow amide rotation: two N-methyls', h: [[8.02, 1], [2.96, 3], [2.88, 3]], c: [[162.62, 1], [36.50, 1], [31.45, 1]] },
  { name: '1,4-dioxane', smiles: 'C1COCCO1', solvent: 'CDCl3', source: 'F2010', h: [[3.71, 8]], c: [[67.14, 4]] },
  { name: 'DME', smiles: 'COCCOC', solvent: 'CDCl3', source: 'F2010', h: [[3.55, 4], [3.40, 6]], c: [[71.84, 2], [59.08, 2]] },
  { name: 'ethanol', smiles: 'CCO', solvent: 'CDCl3', source: 'F2010', tests: 'enantiotopic CH₂ stays one signal', h: [[3.72, 2], [1.25, 3]], c: [[58.28, 1], [18.41, 1]] },
  { name: 'ethyl acetate', smiles: 'CCOC(C)=O', solvent: 'CDCl3', source: 'F2010', h: [[4.12, 2], [2.05, 3], [1.26, 3]], c: [[171.36, 1], [60.49, 1], [21.04, 1], [14.19, 1]] },
  { name: 'ethylene glycol', smiles: 'OCCO', solvent: 'CDCl3', source: 'F2010', h: [[3.76, 4]], c: [[63.79, 2]] },
  { name: 'methanol', smiles: 'CO', solvent: 'CDCl3', source: 'F2010', h: [[3.49, 3]], c: [[50.41, 1]] },
  { name: 'nitromethane', smiles: 'C[N+](=O)[O-]', solvent: 'CDCl3', source: 'F2010', h: [[4.33, 3]], c: [[62.50, 1]] },
  { name: 'n-hexane', smiles: 'CCCCCC', solvent: 'CDCl3', source: 'F2010', tests: 'C₂ / mirror symmetry of a chain', h: [[1.26, 4], [1.26, 4], [0.88, 6]], c: [[31.64, 2], [22.70, 2], [14.14, 2]] },
  { name: 'n-pentane', smiles: 'CCCCC', solvent: 'CDCl3', source: 'F2010', h: [[1.27, 4], [1.27, 2], [0.88, 6]], c: [[34.16, 1], [22.38, 2], [14.08, 2]] },
  { name: '2-propanol', smiles: 'CC(C)O', solvent: 'CDCl3', source: 'F2010', tests: 'enantiotopic methyls stay one signal', h: [[4.04, 1], [1.22, 6]], c: [[64.50, 1], [25.14, 2]] },
  { name: 'pyridine', smiles: 'c1ccncc1', solvent: 'CDCl3', source: 'F2010', h: [[8.62, 2], [7.68, 1], [7.29, 2]], c: [[149.90, 2], [135.96, 1], [123.75, 2]] },
  { name: 'THF', smiles: 'C1CCOC1', solvent: 'CDCl3', source: 'F2010', h: [[3.76, 4], [1.85, 4]], c: [[67.97, 2], [25.62, 2]] },
  { name: 'toluene', smiles: 'Cc1ccccc1', solvent: 'CDCl3', source: 'F2010', h: [[7.25, 2], [7.17, 2], [7.17, 1], [2.36, 3]], c: [[137.89, 1], [129.07, 2], [128.26, 2], [125.33, 1], [21.46, 1]] },
  { name: 'triethylamine', smiles: 'CCN(CC)CC', solvent: 'CDCl3', source: 'F2010', h: [[2.53, 6], [1.03, 9]], c: [[46.25, 3], [11.61, 3]] },
  { name: '18-crown-6', smiles: 'C1COCCOCCOCCOCCOCCO1', solvent: 'CDCl3', source: 'F2010', tests: 'large ring symmetry', h: [[3.67, 24]], c: [[70.55, 12]] },
  { name: 'propane', smiles: 'CCC', solvent: 'CDCl3', source: 'F2010', tests: 'homotopic CH₂', h: [[1.32, 2], [0.90, 6]], c: [[16.63, 2], [16.37, 1]] },
  { name: 'propylene', smiles: 'C=CC', solvent: 'CDCl3', source: 'F2010', tests: '=CH₂ cis and trans H are different', h: [[5.83, 1], [5.03, 1], [4.94, 1], [1.73, 3]], c: [[133.91, 1], [115.74, 1], [19.50, 1]] },
  { name: 'pyrrolidine', smiles: 'C1CCNC1', solvent: 'CDCl3', source: 'F2010', h: [[2.87, 4], [1.68, 4]], c: [[46.93, 2], [25.56, 2]] },
  { name: 'pyrrole', smiles: 'c1cc[nH]c1', solvent: 'CDCl3', source: 'F2010', h: [[6.83, 2], [6.26, 2]], c: [[117.77, 2], [107.98, 2]] },
  { name: 'hexamethylbenzene', smiles: 'Cc1c(C)c(C)c(C)c(C)c1C', solvent: 'CDCl3', source: 'F2010', h: [[2.24, 18]], c: [[132.21, 6], [16.98, 6]] },
  // ---- Stadelmann 2022 (CDCl₃, rigid molecules, diastereotopic H assigned) -------------------------
  { name: 'camphor (29)', smiles: 'CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2', solvent: 'CDCl3', source: 'S2022', tests: 'three different methyls, exo/endo CH₂ protons',
    h: [[2.35, 1], [2.09, 1], [1.95, 1], [1.85, 1], [1.68, 1], [1.41, 1], [1.34, 1], [0.96, 3], [0.91, 3], [0.84, 3]],
    c: [[219.7, 1], [57.7, 1], [46.8, 1], [43.3, 1], [43.1, 1], [29.9, 1], [27.1, 1], [19.8, 1], [19.2, 1], [9.3, 1]] },
  { name: 'fenchone (25)', smiles: 'CC1(C)[C@@H]2CC[C@@](C)(C2)C1=O', solvent: 'CDCl3', source: 'S2022', tests: 'diastereotopic gem-dimethyl, three CH₂ pairs',
    h: [[2.14, 1], [1.79, 1], [1.79, 1], [1.71, 1], [1.56, 1], [1.53, 1], [1.39, 1], [1.14, 3], [1.03, 3], [1.03, 3]],
    c: [[223.4, 1], [54.2, 1], [47.4, 1], [45.4, 1], [41.7, 1], [31.9, 1], [25.0, 1], [23.4, 1], [21.7, 1], [14.6, 1]] },
  { name: 'norcamphor (24)', smiles: 'O=C1C[C@@H]2CC[C@H]1C2', solvent: 'CDCl3', source: 'S2022', tests: 'every CH₂ diastereotopic',
    h: [[2.66, 1], [2.60, 1], [2.05, 1], [1.84, 1], [1.82, 1], [1.80, 1], [1.73, 1], [1.55, 1], [1.53, 1], [1.44, 1]],
    c: [[218.3, 1], [49.9, 1], [45.3, 1], [37.7, 1], [35.3, 1], [27.2, 1], [24.2, 1]] },
  { name: 'norbornene (23)', smiles: 'C1C[C@@H]2C[C@H]1C=C2', solvent: 'CDCl3', source: 'S2022', tests: 'mirror symmetry plus exo/endo and syn/anti H',
    h: [[5.99, 2], [2.85, 2], [1.61, 2], [1.31, 1], [1.08, 1], [0.95, 2]], c: [[135.4, 2], [48.5, 1], [41.7, 2], [24.6, 2]] },
  { name: 'endo-dicyclopentadiene (32)', smiles: 'C1C=C[C@@H]2[C@H]1[C@@H]1C[C@H]2C=C1', solvent: 'CDCl3', source: 'S2022', tests: 'no symmetry: 12 different H',
    h: [[5.98, 1], [5.94, 1], [5.49, 1], [5.49, 1], [3.22, 1], [2.88, 1], [2.78, 1], [2.73, 1], [2.18, 1], [1.68, 1], [1.48, 1], [1.30, 1]],
    c: [[136.0, 1], [132.4, 1], [132.1, 1], [132.0, 1], [54.8, 1], [50.4, 1], [46.2, 1], [45.2, 1], [41.2, 1], [34.7, 1]] },
  { name: 'adamantane (31)', smiles: 'C1C2CC3CC1CC(C2)C3', solvent: 'CDCl3', source: 'S2022', tests: 'Td cage symmetry', h: [[1.88, 4], [1.75, 12]], c: [[37.8, 6], [28.3, 4]] },
  { name: 'isophorone (26)', smiles: 'CC1=CC(=O)CC(C)(C)C1', solvent: 'CDCl3', source: 'S2022', tests: 'gem-dimethyl in a ring without a stereocentre',
    h: [[5.88, 1], [2.20, 2], [2.17, 2], [1.94, 3], [1.04, 6]], c: [[199.9, 1], [160.3, 1], [125.6, 1], [50.8, 1], [45.3, 1], [33.6, 1], [28.3, 2], [24.6, 1]] },
  { name: 'cyclopropanecarbonitrile (8)', smiles: 'N#CC1CC1', solvent: 'CDCl3', source: 'S2022', tests: 'ring CH₂ H cis / trans to CN',
    h: [[1.34, 1], [1.08, 2], [1.01, 2]], c: [[122.2, 1], [7.1, 2], [-3.5, 1]] },
  { name: '3,3-dimethylbut-1-ene (6)', smiles: 'C=CC(C)(C)C', solvent: 'CDCl3', source: 'S2022', tests: '=CH₂ cis and trans H',
    h: [[5.85, 1], [4.92, 1], [4.83, 1], [1.01, 9]], c: [[149.9, 1], [108.9, 1], [33.7, 1], [29.2, 3]] },
  { name: '(E)-crotonaldehyde (4)', smiles: 'C/C=C/C=O', solvent: 'CDCl3', source: 'S2022', h: [[9.50, 1], [6.87, 1], [6.15, 1], [2.03, 3]], c: [[193.9, 1], [153.9, 1], [134.7, 1], [18.7, 1]] },
  { name: 'isobutyronitrile (5)', smiles: 'CC(C)C#N', solvent: 'CDCl3', source: 'S2022', tests: 'enantiotopic methyls', h: [[2.70, 1], [1.33, 6]], c: [[123.8, 1], [20.0, 2], [19.8, 1]] },
  { name: 'propionitrile (3)', smiles: 'CCC#N', solvent: 'CDCl3', source: 'S2022', h: [[2.36, 2], [1.30, 3]], c: [[120.7, 1], [10.9, 1], [10.5, 1]] },
  { name: 'nitroethane (2)', smiles: 'CC[N+](=O)[O-]', solvent: 'CDCl3', source: 'S2022', h: [[4.42, 2], [1.59, 3]], c: [[70.5, 1], [12.3, 1]] },
  { name: 'methyl propiolate (7)', smiles: 'COC(=O)C#C', solvent: 'CDCl3', source: 'S2022', h: [[3.81, 3], [2.88, 1]], c: [[153.1, 1], [74.8, 1], [74.5, 1], [53.0, 1]] },
  { name: 'γ-butyrolactone (9)', smiles: 'O=C1CCCO1', solvent: 'CDCl3', source: 'S2022', h: [[4.35, 2], [2.50, 2], [2.27, 2]], c: [[177.6, 1], [68.5, 1], [27.8, 1], [22.2, 1]] },
  { name: 'dihydrofuran-3(2H)-one (10)', smiles: 'O=C1CCOC1', solvent: 'CDCl3', source: 'S2022', h: [[4.25, 2], [3.87, 2], [2.50, 2]], c: [[215.0, 1], [70.6, 1], [66.9, 1], [37.1, 1]] },
  { name: 'tetrahydropyran-4-one (14)', smiles: 'O=C1CCOCC1', solvent: 'CDCl3', source: 'S2022', h: [[3.98, 4], [2.50, 4]], c: [[206.5, 1], [67.9, 2], [42.9, 2]] },
  { name: 'maleic anhydride (12)', smiles: 'O=C1OC(=O)C=C1', solvent: 'CDCl3', source: 'S2022', h: [[7.03, 2]], c: [[164.1, 2], [136.5, 2]] },
  { name: 'nitrobenzene (15)', smiles: '[O-][N+](=O)c1ccccc1', solvent: 'CDCl3', source: 'S2022', h: [[8.24, 2], [7.70, 1], [7.56, 2]], c: [[148.2, 1], [134.6, 1], [129.3, 2], [123.2, 2]] },
  { name: 'phenol (16)', smiles: 'Oc1ccccc1', solvent: 'CDCl3', source: 'S2022', h: [[7.25, 2], [6.94, 1], [6.83, 2]], c: [] },
  { name: 'anisole (17)', smiles: 'COc1ccccc1', solvent: 'CDCl3', source: 'S2022', h: [[7.30, 2], [6.95, 1], [6.91, 2], [3.81, 3]], c: [[159.6, 1], [129.5, 2], [120.7, 1], [113.9, 2], [55.1, 1]] },
  { name: 'benzaldehyde (18)', smiles: 'O=Cc1ccccc1', solvent: 'CDCl3', source: 'S2022', h: [[10.03, 1], [7.89, 2], [7.64, 1], [7.54, 2]], c: [[192.4, 1], [136.4, 1], [134.5, 1], [129.8, 2], [129.0, 2]] },
  { name: 'benzonitrile (19)', smiles: 'N#Cc1ccccc1', solvent: 'CDCl3', source: 'S2022', h: [[7.67, 2], [7.61, 1], [7.48, 2]], c: [[132.8, 1], [132.2, 2], [129.1, 2], [118.9, 1], [112.5, 1]] },
  { name: '4-formylbenzonitrile (27)', smiles: 'O=Cc1ccc(C#N)cc1', solvent: 'CDCl3', source: 'S2022', h: [[10.10, 1], [8.00, 2], [7.85, 2]], c: [[190.6, 1], [138.8, 1], [132.9, 2], [129.9, 2], [117.7, 1], [117.7, 1]] },
  { name: 'mesitaldehyde (28)', smiles: 'Cc1cc(C)c(C=O)c(C)c1', solvent: 'CDCl3', source: 'S2022', h: [[10.57, 1], [6.90, 2], [2.58, 6], [2.32, 3]], c: [[193.0, 1], [143.8, 1], [141.5, 2], [130.5, 1], [130.0, 2], [21.5, 1], [20.5, 2]] },
  { name: 'DMAP (20)', smiles: 'CN(C)c1ccncc1', solvent: 'CDCl3', source: 'S2022', h: [[8.23, 2], [6.49, 2], [3.00, 6]], c: [[154.2, 1], [149.9, 2], [106.6, 2], [39.0, 2]] },
  { name: 'pyrimidine (21)', smiles: 'c1cncnc1', solvent: 'CDCl3', source: 'S2022', h: [[9.25, 1], [8.76, 2], [7.34, 1]], c: [[159.1, 1], [156.9, 2], [121.6, 1]] },
  { name: 'indole (30)', smiles: 'c1ccc2[nH]ccc2c1', solvent: 'CDCl3', source: 'S2022', h: [[7.65, 1], [7.41, 1], [7.22, 1], [7.20, 1], [7.12, 1], [6.56, 1]],
    c: [[135.8, 1], [127.9, 1], [124.1, 1], [122.0, 1], [120.7, 1], [119.8, 1], [111.0, 1], [102.7, 1]] },
  { name: '1-methylnaphthalene (33)', smiles: 'Cc1cccc2ccccc12', solvent: 'CDCl3', source: 'S2022', h: [[8.00, 1], [7.85, 1], [7.71, 1], [7.52, 1], [7.48, 1], [7.37, 1], [7.32, 1], [2.70, 3]],
    c: [[134.3, 1], [133.5, 1], [132.6, 1], [128.5, 1], [126.6, 1], [126.4, 1], [125.7, 1], [125.6, 1], [125.5, 1], [124.1, 1], [19.4, 1]] },
  { name: 'fluorenone (34)', smiles: 'O=C1c2ccccc2-c2ccccc12', solvent: 'CDCl3', source: 'S2022', tests: 'C₂ᵥ symmetry across two rings',
    h: [[7.66, 2], [7.53, 2], [7.49, 2], [7.30, 2]], c: [[193.9, 1], [144.5, 2], [134.7, 2], [134.2, 2], [129.1, 2], [124.4, 2], [120.3, 2]] },
  { name: 'anthracene (35)', smiles: 'c1ccc2cc3ccccc3cc2c1', solvent: 'CDCl3', source: 'S2022', tests: 'D₂ₕ symmetry', h: [[8.43, 2], [8.01, 4], [7.46, 4]], c: [[131.7, 4], [128.2, 4], [126.2, 2], [125.3, 4]] },
];
