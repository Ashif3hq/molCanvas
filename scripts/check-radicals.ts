// Checks radicals, carbenes and charged species: valence limits, implicit hydrogens, info-bar
// explanations and round trips through SMILES, MOL (V2000/V3000), ChemDraw CDXML and the 3D formats.
// Pass a ChemDraw CDXML file with a radical (e.g. RDKit's chemdraw_template5.cdxml) to also check import.
import { readFileSync } from 'node:fs';
import { Molecule, Resources } from 'openchemlib';
import { analyze, fromMolfile, fromSmiles, toMolfile, toSmiles } from '../src/chem';
import { fromCdxml, toCdxml } from '../src/cdxml';
import { build3D, mol2, molfile3D, pdb, xyz } from '../src/modeling';
import { Atom, Doc, emptyDoc, Radical, addAtom, addBond } from '../src/model';
import { explainAtom, valenceLimit, valenceProblems } from '../src/valence';
import { renderDoc, LIGHT } from '../src/render';

Resources.registerFromNodejs();
let bad = 0;
const fail = (name: string, msg: string) => { bad++; console.error('FAIL', name, msg); };
const eq = (name: string, got: unknown, want: unknown) => { if (JSON.stringify(got) !== JSON.stringify(want)) fail(name, `got ${JSON.stringify(got)}, want ${JSON.stringify(want)}`); };

// --- bond limits ---
const atom = (el: string, charge = 0, radical?: Radical): Atom => ({ id: 1, el, x: 0, y: 0, charge, ...(radical ? { radical } : {}) });
eq('C•', valenceLimit(atom('C', 0, 'doublet')), 3);
eq('triplet C', valenceLimit(atom('C', 0, 'triplet')), 2);
eq('singlet C', valenceLimit(atom('C', 0, 'singlet')), 2);
eq('C•+', valenceLimit(atom('C', 1, 'doublet')), 2);
eq('C•−', valenceLimit(atom('C', -1, 'doublet')), 2);
eq('N•', valenceLimit(atom('N', 0, 'doublet')), 2);
eq('N•+', valenceLimit(atom('N', 1, 'doublet')), 3);
eq('nitrene', valenceLimit(atom('N', 0, 'triplet')), 1);
eq('O•', valenceLimit(atom('O', 0, 'doublet')), 1);
eq('O•−', valenceLimit(atom('O', -1, 'doublet')), 0);
eq('Cl•', valenceLimit(atom('Cl', 0, 'doublet')), 6);
eq('H•', valenceLimit(atom('H', 0, 'doublet')), 0);
eq('H triplet impossible', valenceLimit(atom('H', 0, 'triplet')), -1);

// methane carbon with 4 bonds cannot take a radical; a lone H cannot be a triplet
{
  const d = emptyDoc();
  const c = addAtom(d, 'C', 0, 0);
  for (let i = 0; i < 4; i++) addBond(d, c.id, addAtom(d, 'C', 40 * Math.cos(i), 40 * Math.sin(i)).id, 1);
  c.radical = 'doublet';
  eq('CC4 radical refused', valenceProblems(d).map((p) => p.atom.id), [c.id]);
  const h = emptyDoc();
  addAtom(h, 'H', 0, 0).radical = 'triplet';
  eq('H triplet refused', valenceProblems(h).length, 1);
}

// --- species: SMILES in, implicit H, info bar, SMILES/MOL/CDXML out ---
interface Case { smi: string; idx: number; radical?: Radical; charge: number; h: number; says: string }
const CASES: Record<string, Case> = {
  ethylRadical: { smi: 'C[CH2]', idx: 1, radical: 'doublet', charge: 0, h: 2, says: 'primary carbon radical' },
  tertButylRadical: { smi: 'C[C](C)C', idx: 1, radical: 'doublet', charge: 0, h: 0, says: 'tertiary carbon radical' },
  tempo: { smi: 'CC1(C)CCCC(C)(C)N1[O]', idx: 10, radical: 'doublet', charge: 0, h: 0, says: 'oxyl radical' },
  methylene: { smi: '[CH2]', idx: 0, radical: 'triplet', charge: 0, h: 2, says: 'triplet carbene' },
  radicalCation: { smi: 'C[C+]C', idx: 1, radical: 'doublet', charge: 1, h: 0, says: 'radical cation' },
  isopropylCation: { smi: 'C[CH+]C', idx: 1, charge: 1, h: 1, says: 'secondary carbocation' },
  tButylCation: { smi: 'C[C+](C)C', idx: 1, charge: 1, h: 0, says: 'tertiary carbocation' },
  methoxide: { smi: 'C[O-]', idx: 1, charge: -1, h: 0, says: 'oxyanion' },
  ammonium: { smi: 'C[N+](C)(C)C', idx: 1, charge: 1, h: 0, says: 'ammonium cation' },
  chloride: { smi: '[Cl-]', idx: 0, charge: -1, h: 0, says: 'chloride anion' },
};

const idcode = (d: Doc) => { const m = Molecule.fromSmiles(toSmiles(d)); m.ensureHelperArrays(Molecule.cHelperParities); return m.getIDCode(); };
const radicalsOf = (d: Doc) => d.atoms.map((a) => `${a.el}${a.charge || ''}${a.radical ?? ''}`).filter((s) => /doublet|triplet|singlet|-|\d/.test(s)).sort().join(' ');

for (const [name, c] of Object.entries(CASES)) {
  const doc = fromSmiles(c.smi);
  const a = doc.atoms[c.idx];
  eq(`${name} radical`, a.radical, c.radical);
  eq(`${name} charge`, a.charge, c.charge);
  const an = analyze(doc);
  eq(`${name} implicit H`, an.hCount.get(a.id), c.h);
  eq(`${name} valence ok`, valenceProblems(doc).length, 0);
  const ex = explainAtom(doc, a, an.hCount.get(a.id) ?? 0, false);
  if (!ex.detail.includes(c.says)) fail(name, `info bar says “${ex.detail}”, expected “${c.says}”`);
  if (c.radical && !/unpaired|lone pair/.test(ex.detail)) fail(name, `info bar does not mention the electrons: ${ex.detail}`);

  const want = radicalsOf(doc), ref = idcode(doc);
  const viaSmiles = fromSmiles(toSmiles(doc));
  eq(`${name} SMILES round trip`, radicalsOf(viaSmiles), want);
  for (const v3 of [false, true]) {
    const mf = toMolfile(doc, v3);
    if (c.radical && !(v3 ? /RAD=\d/ : /M {2}RAD/).test(mf)) fail(name, `molfile ${v3 ? 'V3000' : 'V2000'} has no radical`);
    const back = fromMolfile(mf);
    eq(`${name} MOL ${v3 ? 'V3000' : 'V2000'} round trip`, radicalsOf(back), want);
    eq(`${name} MOL structure`, idcode(back), ref);
  }
  const cd = toCdxml(doc);
  if (c.radical && !cd.includes(`Radical="${c.radical[0].toUpperCase()}`)) fail(name, 'CDXML has no Radical attribute');
  const { doc: cdBack, notes } = fromCdxml(cd);
  eq(`${name} CDXML round trip`, radicalsOf(cdBack), want);
  eq(`${name} CDXML structure`, idcode(cdBack), ref);
  if (notes.length) fail(name, `CDXML notes: ${notes.join(' ')}`);

  // drawing: one dot for a doublet, two for a carbene; lone pairs only when asked for
  const dots = (svg: string) => (svg.match(/<circle [^>]*r="1.6"/g) ?? []).length;
  const plain = dots(renderDoc(doc, { theme: LIGHT, colorAtoms: true }));
  eq(`${name} radical dots`, plain, c.radical === 'doublet' ? 1 : c.radical ? 2 : 0);
  if (name === 'methoxide') eq('methoxide lone-pair dots', dots(renderDoc(doc, { theme: LIGHT, colorAtoms: true, showLonePairs: true })), 6);

  // 3D: explicit hydrogens must match the radical, and the 3D molfile must keep it
  const m = build3D(doc);
  const heavyH = m.mol.getAllAtoms() - doc.atoms.length;
  const wantH = [...an.hCount.values()].reduce((s, x) => s + x, 0);
  eq(`${name} 3D hydrogens`, heavyH, wantH);
  if (c.radical && !/M {2}RAD/.test(molfile3D(m))) fail(name, '3D molfile has no radical');
  for (const [fmt, text] of [['MOL2', mol2(m)], ['PDB', pdb(m)], ['XYZ', xyz(m)]]) if (!text.trim()) fail(name, `${fmt} empty`);
}

// --- ChemDraw file with a radical ---
const sample = process.argv[2];
if (sample) {
  const { doc, notes } = fromCdxml(readFileSync(sample, 'utf8'));
  const n = doc.atoms.filter((a) => a.radical).length;
  if (!n) fail('ChemDraw sample', 'no radical imported');
  if (notes.some((t) => /radical|electron/i.test(t))) fail('ChemDraw sample', notes.join(' '));
  console.log(`ChemDraw sample: ${n} radical atom(s), notes: ${notes.join(' ') || 'none'}`);
}

if (bad) { console.error(`${bad} failure(s)`); process.exit(1); }
console.log(`radicals and ions: ${Object.keys(CASES).length} species OK`);
