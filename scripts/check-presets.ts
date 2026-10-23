// Validates every preset: it must parse, and polymers/groups must carry '*' atoms.
import { PRESETS } from '../src/presets';
import { fromSmiles, toSmiles, properties } from '../src/chem';

let bad = 0;
for (const p of PRESETS) {
  try {
    const d = fromSmiles(p.smiles);
    const stars = d.atoms.filter((a) => a.el === '*').length;
    if (p.kind === 'polymer' && stars !== 2) throw new Error(`expected 2 stars, got ${stars}`);
    if (p.kind === 'group' && stars !== 1) throw new Error(`expected 1 star, got ${stars}`);
    if (p.kind === 'molecule' && stars) throw new Error('unexpected star');
    const pr = properties(d);
    if (!pr) throw new Error('no properties');
    if (process.argv.includes('-v')) console.log(p.name.padEnd(34), pr.formula.padEnd(16), pr.mw.toFixed(2), toSmiles(d));
  } catch (e) {
    bad++;
    console.error('FAIL', p.cat, p.name, p.smiles, (e as Error).message);
  }
}
console.log(`${PRESETS.length} presets, ${bad} failures`);
process.exit(bad ? 1 : 0);
