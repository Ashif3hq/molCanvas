// Compares predicted ¹H / ¹³C shifts with experimental values (CDCl₃ unless noted) and checks refusals.
// Each experimental shift is matched to the nearest predicted peak. Additive increments are estimates, so the
// limits are generous: ¹H 0.45 ppm (0.6 for polymers), ¹³C 8 ppm. Run with -v to print every peak.
import { fromSmiles } from '../src/chem';
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
