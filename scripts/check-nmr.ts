// Compares predicted ¹H / ¹³C shifts with experimental values (CDCl₃ unless noted) and checks refusals.
// Each experimental shift is matched to the nearest predicted peak. Additive increments are estimates, so the
// limits are generous: ¹H 0.45 ppm (0.6 for polymers), ¹³C 8 ppm. Run with -v to print every peak.
import { fromSmiles, toOCL } from '../src/chem';
import { equivalence, topicity } from '../src/nmr-symmetry';
import { predictNmr, Nucleus } from '../src/nmr';
import { PRESETS } from '../src/presets';

// [SMILES, ¹H shifts, ¹³C shifts]. Exchangeable O–H / N–H protons are left out.
const CASES: [string, number[], number[]][] = [
  ['CCO', [3.69, 1.22], [58.3, 18.4]], // ethanol
  ['CCOC(C)=O', [4.12, 2.05, 1.26], [171.1, 60.5, 21.0, 14.2]], // ethyl acetate
  ['CC(C)=O', [2.17], [206.7, 30.8]], // acetone
  ['Cc1ccccc1', [7.25, 7.17, 2.35], [137.9, 129.1, 128.3, 125.4, 21.5]], // toluene
  ['C=C(C)C(=O)OC', [6.10, 5.55, 3.75, 1.94], [167.8, 136.3, 125.3, 51.7, 18.3]], // methyl methacrylate
  ['C=Cc1ccccc1', [7.40, 7.32, 7.25, 6.72, 5.74, 5.23], [137.6, 136.9, 128.5, 127.8, 126.2, 113.8]], // styrene
  ['COc1ccccc1', [7.29, 6.93, 6.89, 3.80], [159.7, 129.5, 120.7, 114.0, 55.1]], // anisole
  ['O=Cc1ccccc1', [10.0, 7.87, 7.62, 7.52], [192.4, 136.4, 134.4, 129.7, 129.0]], // benzaldehyde
  ['CCOCC', [3.48, 1.21], [65.9, 15.2]], // diethyl ether
  ['C1CCCCC1', [1.43], [27.0]], // cyclohexane
  ['C1CCOC1', [3.76, 1.85], [68.0, 25.7]], // THF
  ['CCN(CC)CC', [2.53, 1.03], [46.3, 11.6]], // triethylamine
  ['CC(=O)c1ccccc1', [7.96, 7.56, 7.46, 2.61], [198.1, 137.1, 133.1, 128.6, 128.3, 26.6]], // acetophenone
  ['c1ccncc1', [8.62, 7.68, 7.29], [149.9, 135.9, 123.8]], // pyridine
  ['CC(C)O', [4.04, 1.22], [64.3, 25.3]], // 2-propanol
  ['CC#N', [2.10], [116.4, 1.9]], // acetonitrile
  ['CCCCCC', [1.26, 0.88], [31.6, 22.7, 14.1]], // hexane
  ['C=CC(=O)OC', [6.40, 6.13, 5.82, 3.76], [166.6, 130.6, 128.2, 51.5]], // methyl acrylate
  ['C=COC(C)=O', [7.27, 4.88, 4.57, 2.14], [167.9, 141.2, 97.5, 20.6]], // vinyl acetate
  ['ClCCl', [5.30], [53.5]], // dichloromethane
  ['ClC(Cl)Cl', [7.26], [77.2]], // chloroform
  ['CN(C)C=O', [8.02, 2.96, 2.88], [162.6, 36.5, 31.4]], // DMF
  ['[O-][N+](=O)c1ccccc1', [8.23, 7.71, 7.56], [148.3, 134.7, 129.4, 123.5]], // nitrobenzene
  ['CC(O)=O', [2.10], [178.0, 20.8]], // acetic acid
  ['CCc1ccccc1', [7.25, 7.18, 2.65, 1.24], [144.2, 128.4, 127.9, 125.7, 28.9, 15.6]], // ethylbenzene
  ['CCOC(=O)/C=C/c1ccccc1', [7.68, 7.52, 7.38, 6.44, 4.26, 1.34], [167.0, 144.6, 134.5, 130.2, 128.9, 128.0, 118.3, 60.5, 14.3]], // ethyl (E)-cinnamate
  ['CC(C)(C)c1ccccc1', [7.39, 7.30, 7.17, 1.32], [151.1, 128.1, 125.4, 125.2, 34.7, 31.4]], // tert-butylbenzene
  // Polymers (repeat units), literature values
  ['*CC(*)(C)C(=O)OC', [3.60, 1.85, 1.02], [177.8, 54.4, 51.8, 44.8, 17.5]], // PMMA
  ['*CC(*)c1ccccc1', [7.05, 6.60, 1.85, 1.45], [145.5, 127.9, 125.7, 44.0, 40.4]], // polystyrene
  ['*CC*', [1.27], [29.8]], // polyethylene
  ['*CC(*)C', [1.60, 1.25, 0.86], [46.0, 28.5, 21.8]], // polypropylene
  ['*CC(*)Cl', [4.45, 2.15], [57.0, 46.5]], // PVC
  ['*CC(*)OC(C)=O', [4.85, 2.00, 1.80], [170.4, 66.8, 39.3, 21.0]], // PVAc
  ['*CCO*', [3.64], [70.5]], // PEO
  ['*OC(C)C(*)=O', [5.16, 1.58], [169.6, 69.0, 16.6]], // PLA
  ['*[Si](C)(C)O*', [0.07], [1.0]], // PDMS
];
const REFUSED = ['CC(=O)[O-].[Na+]', 'CCP(C)C', '*CCO', '[CH3]'];

const verbose = process.argv.includes('-v');
let bad = 0;
const errs: Record<Nucleus, number[]> = { '1H': [], '13C': [] };
for (const [smi, h, c] of CASES) {
  const poly = smi.includes('*');
  for (const [nuc, exp, lim] of [['1H', h, poly ? 0.6 : 0.45], ['13C', c, 8]] as [Nucleus, number[], number][]) {
    const r = predictNmr(fromSmiles(smi), nuc);
    if (!r.ok) { bad++; console.error('FAIL refused', smi, nuc, r.reason); continue; }
    const pred = r.value.peaks.filter((p) => !p.exchangeable);
    const worst: string[] = [];
    for (const e of exp) {
      const near = pred.reduce((b, p) => (Math.abs(p.shift - e) < Math.abs(b.shift - e) ? p : b));
      const d = near.shift - e;
      errs[nuc].push(Math.abs(d));
      if (Math.abs(d) > lim) { bad++; worst.push(`${e}→${near.shift}`); }
    }
    const list = pred.map((p) => `${p.shift}${p.rough ? '~' : ''} ${p.mult}${nuc === '1H' ? ` ${p.count}H` : ''}`).join(', ');
    if (worst.length) console.error(`FAIL ${nuc.padEnd(3)} ${smi.padEnd(26)} off: ${worst.join(' ')} | ${list}`);
    else if (verbose) console.log(`ok   ${nuc.padEnd(3)} ${smi.padEnd(26)} ${list}`);
  }
}
for (const smi of REFUSED) {
  const r = predictNmr(fromSmiles(smi), '1H');
  if (r.ok) { bad++; console.error('FAIL should refuse', smi); } else if (verbose) console.log('ok   refused', smi, '-', r.reason);
}
// First-order multiplets: trans-alkene doublets (J ≈ 16 Hz) in ethyl cinnamate, triplet/quartet for ethyl.
{
  const r = predictNmr(fromSmiles('CCOC(=O)/C=C/c1ccccc1'), '1H');
  const vinyl = r.ok ? r.value.peaks.filter((p) => p.assignment.includes('alkene')) : [];
  if (vinyl.length !== 2 || vinyl.some((p) => p.mult !== 'd' || !(p.couplings[0]?.j > 14))) { bad++; console.error('FAIL cinnamate vinyl H should be trans doublets', vinyl); }
  const et = r.ok ? r.value.peaks.filter((p) => (p.mult === 'q' || p.mult === 't') && !p.assignment.startsWith('Ar')).map((p) => `${p.mult}${p.count}`).sort().join(',') : '';
  if (et !== 'q2,t3') { bad++; console.error('FAIL ethyl multiplets', et); }
}

// Diastereotopic protons: CH₂ next to a stereocentre or in a rigid ring splits into two 1H signals that
// couple to each other; CH₂ in achiral molecules (ethanol, diethyl ether, cyclohexane) stays one signal.
{
  const diaCH2 = (smi: string) => {
    const r = predictNmr(fromSmiles(smi), '1H');
    return r.ok ? r.value.peaks.filter((p) => p.assignment.includes('diastereotopic (')) : [];
  };
  for (const smi of ['CCO', 'CCOCC', 'C1CCCCC1', 'CCOC(C)=O', 'CC(C)CO']) {
    if (diaCH2(smi).length) { bad++; console.error('FAIL CH₂ should stay equivalent', smi); }
  }
  // 2-butanol C3, deoxypodophyllotoxin (no wedges drawn, and with stereo): H-4, H-11, OCH₂O.
  const expect: [string, number][] = [
    ['CC(O)CC', 2], ['COC(=O)C(N)Cc1ccccc1', 2],
    ['COC1=CC(=CC(=C1OC)OC)C2C3C(CC4=CC5=C(C=C24)OCO5)COC3=O', 6],
    ['COC1=CC(=CC(=C1OC)OC)[C@H]2[C@@H]3[C@H](CC4=CC5=C(C=C24)OCO5)COC3=O', 6],
  ];
  for (const [smi, n] of expect) {
    const ps = diaCH2(smi);
    if (ps.length !== n || ps.some((p) => p.count !== 1 || !p.couplings.length)) { bad++; console.error('FAIL diastereotopic CH₂ not split', smi, ps.map((p) => `${p.shift} ${p.mult}`)); }
  }
  // Deoxypodophyllotoxin lactone OCH₂ (H-11) 4.46 / 3.92 and benzylic H-4 3.08 / 2.75 (CDCl₃): each pair at least 0.2 ppm apart.
  const dpt = diaCH2(expect[2][0]);
  for (const key of ['CH₂–OC(=O)', 'Ar–CH₂']) {
    const pair = dpt.filter((p) => p.assignment.startsWith(key)).map((p) => p.shift);
    if (pair.length !== 2 || Math.abs(pair[0] - pair[1]) < 0.2) { bad++; console.error('FAIL deoxypodophyllotoxin', key, pair); }
  }
  // Valine isopropyl: two methyl doublets.
  const val = predictNmr(fromSmiles('CC(C)C(N)C(=O)O'), '1H');
  if (!val.ok || val.value.peaks.filter((p) => p.assignment === 'CH₃, diastereotopic').length !== 2) { bad++; console.error('FAIL valine methyls should be diastereotopic'); }
}

// Stereo-aware equivalence (nmr-symmetry.ts): topicity of CH₂ protons and methyl pairs, ring and cage symmetry.
{
  const eqOf = (smi: string) => {
    const { mol } = toOCL(fromSmiles(smi));
    return { mol, e: equivalence(mol) };
  };
  const ch2 = (smi: string, atom: number) => {
    const { e } = eqOf(smi);
    const h = e.hCls[atom], p = e.hProper[atom], c = e.hConst[atom];
    return topicity(c[0] === c[1], h[0] === h[1], p[0] === p[1]);
  };
  const cases: [string, number, string, string][] = [
    ['ClCCl', 1, 'homotopic', 'dichloromethane'],
    ['CCC', 1, 'homotopic', 'propane C2'],
    ['CCO', 1, 'enantiotopic', 'ethanol CH₂'],
    ['OCc1ccccc1', 1, 'enantiotopic', 'benzyl alcohol CH₂'],
    ['CC(O)CC', 3, 'diastereotopic', '2-butanol C3 (no wedges)'],
    ['C[C@H](O)CC', 3, 'diastereotopic', '(S)-2-butanol C3'],
    ['O[C@@H](C)C[C@@H](C)O', 3, 'diastereotopic', 'meso-2,4-pentanediol C3'],
    ['O[C@@H](C)C[C@H](C)O', 3, 'homotopic', '(R,R)-2,4-pentanediol C3 (C₂ axis)'],
    ['OC1CCCCC1', 2, 'diastereotopic', 'cyclohexanol C2 (cis / trans to OH)'],
    ['C1CCCCC1', 0, 'homotopic', 'cyclohexane'],
    ['C=CC', 0, 'diastereotopic', 'propene =CH₂'],
    ['C1C2CC3CC1CC(C2)C3', 0, 'homotopic', 'adamantane CH₂'],
    ['C1CC2CCC1C2', 0, 'diastereotopic', 'norbornane C2 exo / endo'],
    ['C1CC2CCC1C2', 6, 'homotopic', 'norbornane C7'],
  ];
  for (const [smi, atom, want, what] of cases) {
    const got = ch2(smi, atom);
    if (got !== want) { bad++; console.error(`FAIL topicity ${what}: ${got}, expected ${want}`); }
  }
  // Methyl pairs: diastereotopic next to a stereocentre (valine, menthol, camphor), enantiotopic in 2-propanol,
  // and the two N-methyls of DMF differ through slow amide rotation.
  const distinct = (smi: string, a: number, b: number) => { const { e } = eqOf(smi); return e.cls[a] !== e.cls[b]; };
  const pairs: [string, number, number, boolean, string][] = [
    ['CC(C)[C@H](N)C(=O)O', 0, 2, true, 'L-valine isopropyl methyls'],
    ['CC(C)C(N)C(=O)O', 0, 2, true, 'valine without wedges'],
    ['CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O', 0, 2, true, 'menthol isopropyl methyls'],
    ['CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2', 0, 2, true, 'camphor C8 / C9'],
    ['CC(C)O', 0, 2, false, '2-propanol methyls'],
    ['CN(C)C=O', 0, 2, true, 'DMF N-methyls'],
    ['CCN(CC)C(C)=O', 0, 4, true, 'N,N-diethylacetamide CH₃'],
    ['CCN(CC)CC', 0, 4, false, 'triethylamine CH₃'],
    ['c1ccc2ccccc2c1', 0, 5, false, 'naphthalene C1 / C5'],
  ];
  for (const [smi, a, b, want, what] of pairs) {
    if (distinct(smi, a, b) !== want) { bad++; console.error(`FAIL equivalence ${what}: expected ${want ? 'different' : 'equivalent'}`); }
  }
}

// Every library preset must either predict without throwing, with finite shifts in a sane window, or be refused.
let covered = 0;
for (const pr of PRESETS) {
  try {
    const h = predictNmr(fromSmiles(pr.smiles), '1H'), c = predictNmr(fromSmiles(pr.smiles), '13C');
    if (h.ok && c.ok) {
      covered++;
      const out = [...h.value.peaks.filter((p) => !Number.isFinite(p.shift) || p.shift < -1 || p.shift > 14), ...c.value.peaks.filter((p) => !Number.isFinite(p.shift) || p.shift < -10 || p.shift > 230)];
      if (out.length) { bad++; console.error('FAIL out of range', pr.name, out.map((p) => p.shift)); }
    }
    if (verbose) console.log('     ', pr.name.padEnd(34), h.ok ? h.value.peaks.map((p) => `${p.shift}${p.rough ? '~' : ''}`).join(' ') : `– ${h.ok ? '' : h.reason}`);
  } catch (e) { bad++; console.error('FAIL threw', pr.name, e); }
}
console.log(`${covered}/${PRESETS.length} presets have predicted spectra.`);
const mae = (a: number[]) => (a.reduce((s, v) => s + v, 0) / a.length).toFixed(2);
console.log(`¹H mean abs. error ${mae(errs['1H'])} ppm over ${errs['1H'].length} shifts; ¹³C ${mae(errs['13C'])} ppm over ${errs['13C'].length}.`);
if (bad) { console.error(`${bad} problem(s)`); process.exit(1); }
console.log('All NMR checks passed.');
