// NMR benchmark: scores the predictor against literature shifts (scripts/nmr-benchmark-data.ts).
//
// Two things are measured for every molecule and nucleus:
//   • Signal count: the number of predicted signals against the number of chemically distinct environments.
//     This is the equivalence test (merged diastereotopic protons give too few, missed symmetry too many).
//   • Shift error: both spectra are expanded to one value per nucleus (a 2H signal counts twice), sorted, and
//     paired in order; the mean absolute difference is the error. No assignment is needed, and merging two
//     environments into one averaged signal is penalised.
//
// The run fails if a molecule is refused, if the signal-count score drops below the floor, or if the mean
// errors rise above the ceilings below (so regressions show up here first). Run with -v for every molecule,
// --json to print the summary as JSON.
import { fromSmiles } from '../src/chem';
import { NmrSolvent, Nucleus, predictNmr } from '../src/nmr';
import { BENCHMARK, SOURCES } from './nmr-benchmark-data';

// Ceilings sit a little above the current scores (¹H 0.12, ¹³C 2.2 ppm); every signal count must be right.
const CEILING: Record<Nucleus, number> = { '1H': 0.15, '13C': 2.6 }; // ppm, mean over all nuclei
const COUNT_FLOOR = 1; // fraction of (molecule, nucleus) pairs with the right number of signals

const verbose = process.argv.includes('-v');
const asJson = process.argv.includes('--json');

interface Score { n: number; sum: number; countOk: number; countAll: number; rows: string[] }
const score: Record<Nucleus, Score> = {
  '1H': { n: 0, sum: 0, countOk: 0, countAll: 0, rows: [] },
  '13C': { n: 0, sum: 0, countOk: 0, countAll: 0, rows: [] },
};
let bad = 0;

const expand = (sig: [number, number][]) => sig.flatMap(([d, k]) => Array<number>(k).fill(d)).sort((a, b) => b - a);

for (const m of BENCHMARK) {
  if (!SOURCES[m.source]) { bad++; console.error(`FAIL ${m.name}: unknown source ${m.source}`); }
  for (const [nuc, exp] of [['1H', m.h], ['13C', m.c]] as [Nucleus, [number, number][]][]) {
    if (!exp.length) continue;
    const r = predictNmr(fromSmiles(m.smiles), nuc, m.solvent as NmrSolvent);
    if (!r.ok) { bad++; console.error(`FAIL ${m.name} ${nuc}: refused (${r.reason})`); continue; }
    const peaks = r.value.peaks.filter((p) => !p.exchangeable);
    const pred = expand(peaks.map((p) => [p.shift, p.count]));
    const want = expand(exp);
    const s = score[nuc];
    if (pred.length !== want.length) {
      bad++;
      console.error(`FAIL ${m.name} ${nuc}: predicted ${pred.length} nuclei, literature ${want.length}`);
      continue;
    }
    let sum = 0;
    for (let k = 0; k < want.length; k++) sum += Math.abs(pred[k] - want[k]);
    s.n += want.length; s.sum += sum;
    s.countAll++;
    const countOk = peaks.length === exp.length;
    if (countOk) s.countOk++;
    const mae = sum / want.length;
    const line = `${countOk ? '  ' : '≠ '}${nuc.padEnd(3)} ${m.name.padEnd(28)} signals ${String(peaks.length).padStart(2)}/${String(exp.length).padEnd(2)} MAE ${mae.toFixed(nuc === '1H' ? 2 : 1).padStart(5)}${m.tests ? `  (${m.tests})` : ''}`;
    s.rows.push(line);
    if (verbose || !countOk) console.log(line);
  }
}

const summary = (['1H', '13C'] as Nucleus[]).map((nuc) => {
  const s = score[nuc];
  return { nucleus: nuc, molecules: s.countAll, nuclei: s.n, meanAbsError: s.n ? s.sum / s.n : NaN, signalCountCorrect: s.countAll ? s.countOk / s.countAll : NaN };
});
if (asJson) console.log(JSON.stringify({ molecules: BENCHMARK.length, summary }, null, 2));
for (const x of summary) {
  console.log(`${x.nucleus === '1H' ? '¹H ' : '¹³C'}: mean abs. error ${x.meanAbsError.toFixed(x.nucleus === '1H' ? 3 : 2)} ppm over ${x.nuclei} nuclei in ${x.molecules} molecules; signal count right for ${Math.round(100 * x.signalCountCorrect)}%.`);
  if (x.meanAbsError > CEILING[x.nucleus]) { bad++; console.error(`FAIL ${x.nucleus} mean error above ${CEILING[x.nucleus]} ppm`); }
  if (x.signalCountCorrect < COUNT_FLOOR) { bad++; console.error(`FAIL ${x.nucleus} signal count right for fewer than ${COUNT_FLOOR * 100}%`); }
}
console.log(`${BENCHMARK.length} molecules from ${new Set(BENCHMARK.map((m) => m.source)).size} sources.`);
if (bad) { console.error(`${bad} problem(s)`); process.exit(1); }
console.log('NMR benchmark passed.');
