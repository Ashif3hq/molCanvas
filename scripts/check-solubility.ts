// Checks the group-contribution solubility code against hand-worked Fedors / Hoftyzer–Van Krevelen sums,
// and that structures outside the tables are refused. Run with -v to list every preset.
import { fromSmiles } from '../src/chem';
import { PRESETS } from '../src/presets';
import { solubility } from '../src/solubility';

// [SMILES, Fedors δ, HVK δD, δP, δH] worked by hand from the published group tables (V = Fedors volume).
// Experimental values (MPa½) are given in comments for orientation only; group contribution is a ±10 % estimate.
const CASES: [string, number, number?, number?, number?][] = [
  ['CCO', 25.73, 15.10, 8.39, 18.32], // ethanol; exp. Hansen 15.8 / 8.8 / 19.4
  ['CC(C)=O', 18.56, 14.52, 9.90, 5.07], // acetone; exp. 15.5 / 10.4 / 7.0
  ['Cc1ccccc1', 18.69, 17.64, 1.05, 0], // toluene; exp. 18.0 / 1.4 / 2.0
  ['C1CCCCC1', 16.51, 16.07, 0, 0], // cyclohexane; exp. 16.8 / 0 / 0.2
  ['Clc1ccccc1', 20.35], // chlorobenzene (sp2 halogen correction); exp. Hildebrand 19.4
  ['*CC(*)c1ccccc1', 21.59], // polystyrene
  ['*CC(*)(C)C(=O)OC', 20.32], // PMMA
  ['*OCCOC(=O)c1ccc(cc1)C(=O)*', 25.40, 21.48, 5.82, 10.77], // PET, ester split across the repeat-unit boundary
  ['*NCCCCCC(=O)*', 25.43], // nylon 6
  ['*CC(*)Cl', 22.57], // PVC
  ['*[Si](C)(C)O*', 15.11], // PDMS (Fedors only: Si has no HVK value)
];
const REFUSED = ['c1ccccc1', 'c1ccncc1', 'c1ccc2ccccc2c1', 'CC(=O)[O-].[Na+]', 'CC(=O)Cl', '*CCO', 'CS(C)=O'];

let bad = 0;
const close = (a: number | undefined, b: number) => a !== undefined && Math.abs(a - b) < 0.02;
for (const [smi, hil, d, p, h] of CASES) {
  const r = solubility(fromSmiles(smi));
  if (!r.ok) { bad++; console.error('FAIL refused', smi, r.reason); continue; }
  const s = r.value;
  const okH = close(s.hildebrand, hil);
  const okHansen = d === undefined || (close(s.hansen?.d, d) && close(s.hansen?.p, p!) && close(s.hansen?.h, h!));
  if (!okH || !okHansen) { bad++; console.error('FAIL', smi, s.hildebrand?.toFixed(2), s.hansen, s.groups); }
  else console.log('ok  ', smi.padEnd(30), `δ ${s.hildebrand!.toFixed(2)}`, s.hansen ? `δD ${s.hansen.d.toFixed(2)} δP ${s.hansen.p.toFixed(2)} δH ${s.hansen.h.toFixed(2)}` : `(HVK: ${s.hvkMissing})`);
}
for (const smi of REFUSED) {
  const r = solubility(fromSmiles(smi));
  if (r.ok) { bad++; console.error('FAIL should refuse', smi, r.value.groups); }
  else console.log('ok  ', smi.padEnd(30), 'refused:', r.reason);
}
let covered = 0;
for (const pr of PRESETS) {
  const r = solubility(fromSmiles(pr.smiles));
  if (r.ok) covered++;
  if (process.argv.includes('-v')) console.log(pr.name.padEnd(34), r.ok ? `δ ${r.value.hildebrand!.toFixed(1)}${r.value.hansen ? ` (${r.value.hansen.d.toFixed(1)}, ${r.value.hansen.p.toFixed(1)}, ${r.value.hansen.h.toFixed(1)})` : ''}` : `– ${r.reason}`);
}
console.log(`${CASES.length + REFUSED.length} checks, ${bad} failures; ${covered}/${PRESETS.length} presets have solubility parameters`);
process.exit(bad ? 1 : 0);
