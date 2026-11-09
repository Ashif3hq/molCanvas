// Reference NMR data for the spectrum dialog: functional-group ranges, residual solvent signals, common
// impurities and literature spectra of common polymers. Values are ppm; a pair [lo, hi] is a range.
//
// Residual solvent and impurity shifts (CDCl₃): G.R. Fulmer et al., Organometallics 29 (2010) 2176.
// Group ranges: Pretsch, Bühlmann & Badertscher, Structure Determination of Organic Compounds (2009);
//   Silverstein, Webster & Kiemle, Spectrometric Identification of Organic Compounds, 7th ed. (2005).
// Polymer spectra: typical literature values (Bovey & Mirau, NMR of Polymers, Academic Press, 1996, and
//   the polymer NMR literature); exact positions move with solvent, temperature and tacticity.
import { NmrSolvent } from './nmr';

export type RefValue = number | [number, number];
export type RefCategory = 'group' | 'solvent' | 'impurity' | 'polymer';

export interface NmrRef {
  id: string;
  name: string;
  cat: RefCategory;
  h?: RefValue[];
  c?: RefValue[];
  note?: string;
  /** SMILES of the compound or repeat unit, so a drawn structure can be matched to its literature spectrum. */
  smiles?: string;
}

export const SOLVENT_PEAKS: Record<NmrSolvent, { h: number[]; c: number[]; water: number; label: string }> = {
  CDCl3: { h: [7.26], c: [77.16], water: 1.56, label: 'CDCl₃' },
  'DMSO-d6': { h: [2.50], c: [39.52], water: 3.33, label: 'DMSO-d₆' },
  CD3OD: { h: [3.31], c: [49.00], water: 4.87, label: 'CD₃OD' },
  D2O: { h: [4.79], c: [], water: 4.79, label: 'D₂O' },
  'acetone-d6': { h: [2.05], c: [29.84, 206.26], water: 2.84, label: 'acetone-d₆' },
  C6D6: { h: [7.16], c: [128.06], water: 0.40, label: 'C₆D₆' },
  CD3CN: { h: [1.94], c: [1.32, 118.26], water: 2.13, label: 'CD₃CN' },
  CD2Cl2: { h: [5.32], c: [53.84], water: 1.52, label: 'CD₂Cl₂' },
};

const S = (id: string, name: string, h: RefValue[], c: RefValue[], smiles?: string, note?: string): NmrRef =>
  ({ id, name, cat: 'impurity', h, c, smiles, note: note ?? 'in CDCl₃ · Fulmer 2010' });
const G = (id: string, name: string, h?: RefValue[], c?: RefValue[], note?: string): NmrRef => ({ id, name, cat: 'group', h, c, note });
const P = (id: string, name: string, smiles: string, h: RefValue[], c: RefValue[], note: string): NmrRef => ({ id, name, cat: 'polymer', smiles, h, c, note });

export const NMR_REFS: NmrRef[] = [
  // Functional-group ranges (¹H and ¹³C)
  G('g-ch3', 'Alkyl CH₃', [[0.8, 1.0]], [[10, 25]], 'CH₃ of alkyl chains'),
  G('g-ch2', 'Alkyl CH₂ / CH', [[1.2, 1.6]], [[20, 45]], 'chain CH₂, ring CH₂'),
  G('g-allyl', 'Allylic C–H', [[1.6, 2.2]], [[18, 35]], 'C=C–CH'),
  G('g-acyl', 'α to C=O (ketone, ester, acid)', [[2.0, 2.6]], [[20, 45]], 'CH₃C(=O), CH₂C(=O)'),
  G('g-benzyl', 'Benzylic C–H', [[2.2, 2.8]], [[20, 40]], 'Ar–CH₃ 2.3, Ar–CH₂ 2.6'),
  G('g-alkyneh', 'Terminal alkyne ≡C–H', [[1.8, 3.1]], [[65, 90]], '¹³C: ≡CH 68–78, ≡C–R 80–90'),
  G('g-nch', 'Amine N–CH', [[2.2, 2.9]], [[30, 60]], 'N–CH₃ 2.2–2.4'),
  G('g-amidench', 'Amide N–CH', [[2.7, 3.5]], [[26, 45]], 'C(=O)N–CH₃ 2.8–3.0'),
  G('g-och', 'Ether / alcohol O–CH', [[3.3, 4.0]], [[50, 80]], 'O–CH₃ 3.3; O–CH₂ 3.4–3.7'),
  G('g-ome-ester', 'Methyl ester O–CH₃', [[3.6, 3.9]], [[51, 53]], 'C(=O)O–CH₃'),
  G('g-oester', 'Ester O–CH (alkyl side)', [[3.9, 5.2]], [[60, 75]], 'C(=O)O–CH₂ 4.1, C(=O)O–CH 4.9–5.1'),
  G('g-arome', 'Aryl ether O–CH₃', [[3.7, 3.9]], [[55, 56]], 'anisole-type'),
  G('g-halide', 'C–Cl / C–Br / C–I', [[3.1, 3.7]], [[-5, 50]], 'CH₂Cl 3.5, CH₂Br 3.4, CH₂I 3.2'),
  G('g-alkene', 'Alkene =C–H', [[4.5, 6.5]], [[100, 150]], 'vinyl, acrylate, enone'),
  G('g-vinylo', 'Enol ether / vinyl ester =CH–O', [[6.2, 7.3]], [[140, 155]], 'vinyl acetate 7.27'),
  G('g-acrylate', 'Acrylate / methacrylate =CH₂', [[5.5, 6.5]], [[125, 137]], 'residual monomer check'),
  G('g-arom', 'Aromatic C–H', [[6.5, 8.5]], [[110, 160]], 'benzene 7.36'),
  G('g-hetarom', 'Heteroaromatic C–H', [[6.2, 9.3]], [[100, 155]], 'pyridine α 8.6, furan β 6.4'),
  G('g-cho', 'Aldehyde CHO', [[9.4, 10.1]], [[190, 205]]),
  G('g-cooh', 'Carboxylic acid O–H', [[10, 13]], [[170, 185]], 'broad, exchangeable'),
  G('g-oh', 'Alcohol O–H', [[1, 5]], undefined, 'broad; position depends on concentration and water'),
  G('g-phoh', 'Phenol O–H', [[4, 8]], undefined, 'broad; 8.5–10 in DMSO-d₆'),
  G('g-nh', 'Amine N–H', [[0.5, 3]], undefined, 'broad, exchangeable'),
  G('g-amidenh', 'Amide N–H', [[5, 8.5]], undefined, 'broad; 7.5–9 in DMSO-d₆'),
  G('g-cn', 'Nitrile C≡N', undefined, [[115, 125]]),
  G('g-ester', 'Ester / lactone C=O', undefined, [[165, 178]], 'aryl esters 165–167'),
  G('g-amide', 'Amide / lactam C=O', undefined, [[165, 175]]),
  G('g-acid', 'Carboxylic acid C=O', undefined, [[170, 185]]),
  G('g-ketone', 'Ketone C=O', undefined, [[195, 220]], 'aryl ketones 195–200'),
  G('g-carbonate', 'Carbonate / carbamate / urea C=O', undefined, [[150, 160]]),
  G('g-cf3', 'CF₃ (quartet in ¹³C)', undefined, [[120, 128]], '¹J(C,F) ≈ 270–290 Hz'),
  G('g-sime', 'Si–CH₃', [[-0.1, 0.3]], [[-3, 3]], 'silicones, TMS ethers'),

  // Residual solvent signals and common impurities
  { id: 's-cdcl3', name: 'CHCl₃ in CDCl₃', cat: 'solvent', h: [7.26], c: [77.16], note: 'residual solvent · Fulmer 2010' },
  { id: 's-dmso', name: 'DMSO-d₅ in DMSO-d₆', cat: 'solvent', h: [2.50], c: [39.52], note: 'residual solvent; water 3.33' },
  { id: 's-cd3od', name: 'CHD₂OD in CD₃OD', cat: 'solvent', h: [3.31], c: [49.00], note: 'residual solvent; water 4.87' },
  { id: 's-d2o', name: 'HDO in D₂O', cat: 'solvent', h: [4.79], note: 'residual solvent' },
  { id: 's-acetone', name: 'Acetone-d₅ in acetone-d₆', cat: 'solvent', h: [2.05], c: [29.84, 206.26], note: 'residual solvent; water 2.84' },
  { id: 's-c6d6', name: 'C₆D₅H in C₆D₆', cat: 'solvent', h: [7.16], c: [128.06], note: 'residual solvent; water 0.40' },
  { id: 's-cd3cn', name: 'CHD₂CN in CD₃CN', cat: 'solvent', h: [1.94], c: [1.32, 118.26], note: 'residual solvent; water 2.13' },
  { id: 's-cd2cl2', name: 'CHDCl₂ in CD₂Cl₂', cat: 'solvent', h: [5.32], c: [53.84], note: 'residual solvent; water 1.52' },
  { id: 's-water', name: 'Water in CDCl₃', cat: 'solvent', h: [1.56], note: 'Fulmer 2010' },
  { id: 's-tms', name: 'TMS (reference)', cat: 'solvent', h: [0.00], c: [0.00], note: 'δ = 0 by definition' },
  S('i-acetone', 'Acetone', [2.17], [207.07, 30.92], 'CC(C)=O'),
  S('i-mecn', 'Acetonitrile', [2.10], [116.43, 1.89], 'CC#N'),
  S('i-acoh', 'Acetic acid', [2.10], [175.99, 20.81], 'CC(O)=O'),
  S('i-benzene', 'Benzene', [7.36], [128.37], 'c1ccccc1'),
  S('i-tbuoh', 'tert-Butanol', [1.28], [69.15, 31.25], 'CC(C)(C)O'),
  S('i-chcl3', 'Chloroform', [7.26], [77.36], 'ClC(Cl)Cl'),
  S('i-chex', 'Cyclohexane', [1.43], [26.94], 'C1CCCCC1'),
  S('i-dcm', 'Dichloromethane', [5.30], [53.52], 'ClCCl'),
  S('i-et2o', 'Diethyl ether', [3.48, 1.21], [65.91, 15.20], 'CCOCC'),
  S('i-dmf', 'DMF', [8.02, 2.96, 2.88], [162.62, 36.50, 31.39], 'CN(C)C=O'),
  S('i-dmso', 'DMSO', [2.62], [40.76], 'CS(C)=O'),
  S('i-etoh', 'Ethanol', [3.72, 1.25], [58.28, 18.41], 'CCO'),
  S('i-etoac', 'Ethyl acetate', [4.12, 2.05, 1.26], [171.36, 60.49, 21.04, 14.19], 'CCOC(C)=O'),
  S('i-grease', 'Grease (hydrocarbon)', [1.26, 0.86], [29.76]),
  S('i-hexane', 'n-Hexane', [1.26, 0.88], [31.64, 22.70, 14.14], 'CCCCCC'),
  S('i-meoh', 'Methanol', [3.49], [50.41], 'CO'),
  S('i-ipa', '2-Propanol', [4.04, 1.22], [64.50, 25.14], 'CC(C)O'),
  S('i-pyr', 'Pyridine', [8.62, 7.68, 7.29], [149.90, 135.96, 123.75], 'c1ccncc1'),
  S('i-silicone', 'Silicone grease', [0.07], [1.04]),
  S('i-thf', 'THF', [3.76, 1.85], [67.97, 25.62], 'C1CCOC1'),
  S('i-toluene', 'Toluene', [7.25, 7.17, 2.36], [137.89, 129.07, 128.26, 125.33, 21.46], 'Cc1ccccc1'),
  S('i-et3n', 'Triethylamine', [2.53, 1.03], [46.25, 11.61], 'CCN(CC)CC'),

  // Polymers
  P('p-pe', 'Polyethylene (PE)', '*CC*', [1.30], [29.9], 'TCE-d₂ or ODCB-d₄ at 120 °C; CH₂ only'),
  P('p-pp', 'Polypropylene (PP)', '*CC(*)C', [1.60, 1.25, 0.86], [46.0, 28.5, 21.8], 'isotactic; CH₂ protons split 1.25/0.85, CH₃ mmmm 21.8'),
  P('p-ps', 'Polystyrene (PS)', '*CC(*)c1ccccc1', [7.05, 6.60, 1.85, 1.45], [145.5, 127.9, 125.7, [40, 46], 40.4], 'atactic, CDCl₃; ortho-H 6.3–6.7 shielded by neighbouring rings'),
  P('p-pmma', 'Poly(methyl methacrylate) (PMMA)', '*CC(*)(C)C(=O)OC', [3.60, [1.8, 1.9], 1.21, 1.02, 0.84], [177.8, 54.4, 51.8, 44.8, 18.8, 16.5], 'CDCl₃; α-CH₃ triads mm 1.21, mr 1.02, rr 0.84'),
  P('p-pvc', 'Poly(vinyl chloride) (PVC)', '*CC(*)Cl', [[4.3, 4.6], [2.0, 2.3]], [[56, 58], [45, 47]], 'THF-d₈ or ODCB-d₄; split by tacticity'),
  P('p-pvac', 'Poly(vinyl acetate) (PVAc)', '*CC(*)OC(C)=O', [4.85, 2.00, 1.75], [170.4, 66.8, 39.3, 21.0], 'CDCl₃'),
  P('p-pva', 'Poly(vinyl alcohol) (PVA)', '*CC(*)O', [[4.2, 4.7], 3.85, [1.3, 1.5]], [[64, 69], [44, 46]], 'DMSO-d₆; OH triads 4.2–4.7'),
  P('p-peo', 'Poly(ethylene oxide) / PEG', '*CCO*', [3.64], [70.5], 'CDCl₃'),
  P('p-pla', 'Poly(lactic acid) (PLA)', '*OC(C)C(*)=O', [5.16, 1.58], [169.6, 69.0, 16.6], 'CDCl₃; CH quartet, CH₃ doublet'),
  P('p-pcl', 'Poly(ε-caprolactone) (PCL)', '*OCCCCCC(*)=O', [4.06, 2.30, 1.65, 1.38], [173.5, 64.1, 34.1, 28.3, 25.5, 24.6], 'CDCl₃'),
  P('p-pdms', 'Poly(dimethylsiloxane) (PDMS)', '*[Si](C)(C)O*', [0.07], [1.0], 'CDCl₃'),
  P('p-pi', 'cis-1,4-Polyisoprene (natural rubber)', '*CC(C)=CC*', [5.12, 2.04, 1.68], [135.2, 125.0, 32.2, 26.4, 23.4], 'CDCl₃'),
  P('p-pb', '1,4-Polybutadiene', '*CC=CC*', [5.40, 2.05], [130.0, 129.5, 32.7, 27.4], 'CDCl₃; trans CH₂ 32.7, cis 27.4; 1,2-units =CH₂ 4.9–5.0'),
  P('p-pet', 'Poly(ethylene terephthalate) (PET)', '*OCCOC(=O)c1ccc(cc1)C(*)=O', [8.10, 4.75], [165.5, 133.8, 129.8, 63.5], 'CDCl₃/TFA-d'),
  P('p-pc', 'Bisphenol A polycarbonate (PC)', '*Oc1ccc(cc1)C(C)(C)c1ccc(cc1)OC(*)=O', [7.25, 7.16, 1.69], [152.0, 149.0, 148.5, 128.0, 120.3, 42.6, 30.9], 'CDCl₃'),
  P('p-pan', 'Polyacrylonitrile (PAN)', '*CC(*)C#N', [[3.0, 3.2], [2.0, 2.1]], [120.0, 32.5, 27.0], 'DMSO-d₆'),
  P('p-pvdf', 'Poly(vinylidene fluoride) (PVDF)', '*CC(*)(F)F', [2.90, 2.30], [120.0, 43.5], 'DMF-d₇ or acetone-d₆; 2.3 = head-to-head defects; CF₂ is a triplet in ¹³C'),
];

export const REF_CATEGORIES: { id: RefCategory | 'all'; label: string }[] = [
  { id: 'all', label: 'All' },
  { id: 'group', label: 'Functional groups' },
  { id: 'solvent', label: 'Solvents' },
  { id: 'impurity', label: 'Impurities' },
  { id: 'polymer', label: 'Polymers' },
];

export const fmtRef = (vs: RefValue[] | undefined, d: number) =>
  (vs ?? []).map((v) => (Array.isArray(v) ? `${v[0].toFixed(d === 2 ? 1 : 0)}–${v[1].toFixed(d === 2 ? 1 : 0)}` : v.toFixed(d))).join(', ');
