// Checks in-chain repeat-unit brackets (Bracket tool): a block copolymer *–[styrene]ₙ–[MMA]ₘ–*
// gets one bracket pair per block with its own count, each block's M₀, and a CDXML round trip.
import { Resources } from 'openchemlib';
import { fromSmiles, properties } from '../src/chem';
import { fromCdxml, toCdxml } from '../src/cdxml';
import { addAtom, addBond, bracketProblem, Doc, repeatUnitEnds, subDoc, validBrackets } from '../src/model';
import { bracketMarks, LIGHT, renderDoc } from '../src/render';

Resources.registerFromNodejs();
let bad = 0;
const fail = (name: string, msg: string) => { bad++; console.error('FAIL', name, msg); };
const eq = (name: string, got: unknown, want: unknown) => { if (JSON.stringify(got) !== JSON.stringify(want)) fail(name, `got ${JSON.stringify(got)}, want ${JSON.stringify(want)}`); };
const near = (name: string, got: number, want: number, tol = 0.05) => { if (Math.abs(got - want) > tol) fail(name, `got ${got}, want ${want}`); };

// *–CH2–CH(Ph)–CH2–C(CH3)(CO2CH3)–*  (SMILES atom order is kept)
const doc = fromSmiles('*CC(c1ccccc1)CC(C)(C(=O)OC)*');
const ids = doc.atoms.map((a) => a.id);
const styrene = new Set(ids.slice(1, 9)); // CH2, CH and the phenyl ring
const mma = new Set(ids.slice(9, 16)); // CH2, C, CH3, C(=O)OCH3

eq('styrene block ok', bracketProblem(doc, styrene), null);
eq('mma block ok', bracketProblem(doc, mma), null);
eq('whole chain is not a block', bracketProblem(doc, new Set([...styrene, ...mma, ids[0], ids[16]])) !== null, true);
eq('half a ring is refused', bracketProblem(doc, new Set(ids.slice(1, 5))) !== null, true);
eq('two * ends: one automatic unit before brackets', repeatUnitEnds(doc).size, 2);

doc.brackets = [
  { id: doc.nextId++, atoms: [...styrene], sub: 'n' },
  { id: doc.nextId++, atoms: [...mma], sub: 'm' },
];
eq('two valid blocks', validBrackets(doc).length, 2);
eq('* ends are plain asterisks once blocks are bracketed', repeatUnitEnds(doc).size, 0);
const marks = bracketMarks(doc);
eq('four brackets', marks.length, 4);
eq('one count per block', marks.filter((m) => m.sub).map((m) => m.sub!.text).sort(), ['m', 'n']);
const shared = marks.filter((m) => marks.filter((x) => x.bond === m.bond).length === 2);
eq('blocks meeting at one bond get two brackets ][', shared.length, 2);
if (shared.length === 2) {
  const gap = Math.hypot(shared[0].c.x - shared[1].c.x, shared[0].c.y - shared[1].c.y);
  if (gap < 8) fail('][ gap', `brackets on the shared bond overlap (${gap.toFixed(1)} px apart)`);
}
const svg = renderDoc(doc, { theme: LIGHT, colorAtoms: true });
eq('rendered counts', [...svg.matchAll(/data-sub="\d+">([^<]+)</g)].map((m) => m[1]).sort(), ['m', 'n']);

// each block's repeat-unit mass: styrene C8H8 104.15, MMA C5H8O2 100.12
const unit = (d: Doc, atoms: Set<number>) => {
  const u = subDoc(d, atoms);
  for (const c of validBrackets(d).find((b) => b.group.atoms.every((id) => atoms.has(id)))!.crossings) addBond(u, c.inner.id, addAtom(u, '*', c.outer.x, c.outer.y).id);
  return properties(u)!;
};
const ps = unit(doc, styrene), pm = unit(doc, mma);
eq('styrene formula', ps.formula, 'C8H8');
near('styrene M0', ps.mw, 104.15);
eq('MMA formula', pm.formula, 'C5H8O2');
near('MMA M0', pm.mw, 100.12);
eq('whole drawing flagged as bracketed', properties(doc)!.bracketBlocks, 2);

// CDXML round trip keeps both blocks and their counts
const xml = toCdxml(doc);
eq('two SRU groups written', (xml.match(/<bracketedgroup /g) ?? []).length, 2);
const back = fromCdxml(xml).doc;
const vb = validBrackets(back);
eq('two blocks read back', vb.length, 2);
eq('counts read back', vb.map((b) => b.group.sub).sort(), ['m', 'n']);
eq('block sizes read back', vb.map((b) => b.group.atoms.length).sort(), [7, 8]);

// a homopolymer drawn with two * ends still gets the automatic bracket and n
const hp = fromSmiles('*CC(c1ccccc1)*');
eq('homopolymer automatic brackets', bracketMarks(hp).length, 2);

if (bad) { console.error(`${bad} bracket check(s) failed`); process.exit(1); }
console.log('brackets ok');
