// Checks MOL (V2000/V3000) round trips and the 3D modelling formats (MOL, SDF, MOL2, PDB, XYZ).
// Pass an output directory to also write every file there, e.g. for validation with RDKit or Open Babel.
import { mkdirSync, writeFileSync } from 'node:fs';
import { Molecule, Resources } from 'openchemlib';
import { fromMolfile, fromSmiles, toMolfile, toOCL, toSmiles } from '../src/chem';
import { build3D, mol2, molfile3D, pdb, sdf3D, xyz } from '../src/modeling';

Resources.registerFromNodejs();
const outDir = process.argv[2];
if (outDir) mkdirSync(outDir, { recursive: true });

const CASES: Record<string, string> = {
  alanine: 'C[C@H](N)C(=O)O',
  ibuprofen: 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O',
  nitrobenzene: 'c1ccccc1[N+](=O)[O-]',
  acetate: 'CC(=O)[O-]',
  nmethylacetamide: 'CC(=O)NC',
  tetramethylammonium: 'C[N+](C)(C)C',
  transStilbene: 'C(/c1ccccc1)=C\\c1ccccc1',
  pmma: '*CC(*)(C)C(=O)OC',
  pet: '*OCCOC(=O)c1ccc(cc1)C(=O)*',
};

let bad = 0;
const fail = (name: string, msg: string) => { bad++; console.error('FAIL', name, msg); };
const idcode = (m: InstanceType<typeof Molecule>) => { m.ensureHelperArrays(Molecule.cHelperParities); return m.getIDCode(); };

for (const [name, smi] of Object.entries(CASES)) {
  const doc = fromSmiles(smi);
  const ref = toSmiles(doc);
  // 2D molfiles, both versions, must read back to the same structure (including * ends and stereo).
  for (const v3 of [false, true]) {
    const text = toMolfile(doc, v3);
    if (text.includes(' ? ')) fail(name, 'molfile still contains ? atoms');
    const back = toSmiles(fromMolfile(text));
    if (back !== ref) fail(name, `${v3 ? 'V3000' : 'V2000'} round trip gave ${back}, expected ${ref}`);
  }
  // 3D: the conformer must keep constitution and stereo (stars capped with H).
  const capped = toOCL(doc).mol;
  for (let i = 0; i < capped.getAllAtoms(); i++) if (capped.getAtomicNo(i) === 0) capped.setAtomicNo(i, 1);
  const m = build3D(doc);
  const n = m.mol.getAllAtoms();
  for (const v3 of [false, true]) {
    if (idcode(Molecule.fromMolfile(molfile3D(m, v3))) !== idcode(capped)) fail(name, `3D ${v3 ? 'V3000' : 'V2000'} changed the structure or stereo`);
  }
  const sdf = sdf3D(m);
  if (!sdf.trimEnd().endsWith('$$$$') || !sdf.includes('> <SMILES>')) fail(name, 'SDF record malformed');
  const m2 = mol2(m);
  const atomLines = m2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>BOND')[0].trim().split('\n');
  if (atomLines.length !== n) fail(name, `MOL2 has ${atomLines.length} atoms, expected ${n}`);
  const hetatm = pdb(m).split('\n').filter((l) => l.startsWith('HETATM'));
  if (hetatm.length !== n || hetatm.some((l) => l.length !== 80)) fail(name, 'PDB HETATM records malformed');
  if (xyz(m).trim().split('\n').length !== n + 2) fail(name, 'XYZ malformed');
  // Bond lengths after MMFF94s+ must be physical (0.9–2.0 Å).
  for (let b = 0; b < m.mol.getAllBonds(); b++) {
    const a = m.mol.getBondAtom(0, b), c = m.mol.getBondAtom(1, b);
    const d = Math.hypot(m.mol.getAtomX(a) - m.mol.getAtomX(c), m.mol.getAtomY(a) - m.mol.getAtomY(c), m.mol.getAtomZ(a) - m.mol.getAtomZ(c));
    if (d < 0.9 || d > 2.0) { fail(name, `bond ${b} is ${d.toFixed(2)} Å`); break; }
  }
  if (!m.minimised) fail(name, 'MMFF94s+ minimisation did not converge');
  console.log('ok  ', name.padEnd(20), `${n} atoms`, m.energy !== undefined ? `E(MMFF94s+) ${m.energy.toFixed(2)} kcal/mol` : '', m.capped ? `(${m.capped} * capped with H)` : '');
  if (outDir) {
    const w = (ext: string, text: string) => writeFileSync(`${outDir}/${name}.${ext}`, text);
    w('smi', `${smi}\n`); w('2d.mol', toMolfile(doc)); w('2d.v3.mol', toMolfile(doc, true));
    w('mol', molfile3D(m)); w('v3.mol', molfile3D(m, true)); w('sdf', sdf); w('mol2', m2); w('pdb', pdb(m)); w('xyz', xyz(m));
  }
}
console.log(`${Object.keys(CASES).length} structures, ${bad} failures`);
process.exit(bad ? 1 : 0);
