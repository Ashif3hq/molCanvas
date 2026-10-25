// Chemistry layer: converts between the MolCanvas document and OpenChemLib molecules
// for valence, rings, formulas, properties, SMILES and MOL files.
import { Molecule, MoleculeProperties } from 'openchemlib';
import { BOND, Doc, emptyDoc, fragments, Point, repeatUnitEnds } from './model';
import { valenceProblems } from './valence';

type Mol = InstanceType<typeof Molecule>;

const STAR_PLACEHOLDER = 117; // Ts: never occurs in real structures, used to carry '*' through the SMILES parser

export function atomicNo(el: string): number {
  if (el === '*') return 0;
  const n = Molecule.getAtomicNoFromLabel(el);
  return n > 0 && n < 119 ? n : -1;
}

export function isElement(el: string): boolean {
  return atomicNo(el) > 0;
}

export interface OclBuild {
  mol: Mol;
  index: Map<number, number>; // atom id -> OCL atom index
  ids: number[]; // OCL atom index -> atom id
  bondIndex: Map<number, number>; // bond id -> OCL bond index
}

export function toOCL(doc: Doc): OclBuild {
  const mol = new Molecule(doc.atoms.length + 2, doc.bonds.length + 2);
  const index = new Map<number, number>();
  const ids: number[] = [];
  const bondIndex = new Map<number, number>();
  for (const a of doc.atoms) {
    const z = atomicNo(a.el);
    const i = mol.addAtom(z < 0 ? 6 : z);
    mol.setAtomX(i, a.x / BOND);
    mol.setAtomY(i, a.y / BOND);
    if (a.charge) mol.setAtomCharge(i, a.charge);
    if (a.iso) mol.setAtomMass(i, a.iso);
    index.set(a.id, i);
    ids.push(a.id);
  }
  for (const b of doc.bonds) {
    const ia = index.get(b.a), ib = index.get(b.b);
    if (ia === undefined || ib === undefined) continue;
    let type = b.order === 3 ? Molecule.cBondTypeTriple : b.order === 2 ? Molecule.cBondTypeDouble : Molecule.cBondTypeSingle;
    if (b.order === 1 && b.stereo === 'wedge') type = Molecule.cBondTypeUp;
    if (b.order === 1 && b.stereo === 'hash') type = Molecule.cBondTypeDown;
    const bi = mol.addBond(ia, ib);
    mol.setBondType(bi, type);
    bondIndex.set(b.id, bi);
  }
  mol.ensureHelperArrays(Molecule.cHelperRings);
  return { mol, index, ids, bondIndex };
}

export interface Analysis {
  hCount: Map<number, number>;
  badValence: Set<number>;
  ringCenter: Map<number, Point>; // bond id -> centre of the ring used to place its second line
}

export function analyze(doc: Doc): Analysis {
  const hCount = new Map<number, number>();
  const badValence = new Set<number>();
  const ringCenter = new Map<number, Point>();
  if (!doc.atoms.length) return { hCount, badValence, ringCenter };
  const { mol, ids, bondIndex } = toOCL(doc);
  for (let i = 0; i < ids.length; i++) {
    const z = mol.getAtomicNo(i);
    if (z === 0) continue;
    hCount.set(ids[i], mol.getImplicitHydrogens(i));
  }
  for (const p of valenceProblems(doc)) badValence.add(p.atom.id);
  const rings = mol.getRingSet();
  const byAtom = new Map<number, Point>();
  for (const a of doc.atoms) byAtom.set(a.id, a);
  // For each bond, prefer the smallest ring, and among equals the one with most double bonds.
  const ringInfo: { bonds: number[]; center: Point; score: number }[] = [];
  for (let r = 0; r < rings.getSize(); r++) {
    const atoms = rings.getRingAtoms(r);
    const pts = atoms.map((i) => byAtom.get(ids[i])!);
    const center = { x: pts.reduce((s, p) => s + p.x, 0) / pts.length, y: pts.reduce((s, p) => s + p.y, 0) / pts.length };
    const rb = rings.getRingBonds(r);
    let doubles = 0;
    for (const bi of rb) if (mol.getBondOrder(bi) === 2) doubles++;
    ringInfo.push({ bonds: rb, center, score: rings.getRingSize(r) * 10 - doubles });
  }
  ringInfo.sort((a, b) => a.score - b.score);
  const oclToBond = new Map<number, number>();
  for (const [id, bi] of bondIndex) oclToBond.set(bi, id);
  for (const ri of ringInfo) {
    for (const bi of ri.bonds) {
      const id = oclToBond.get(bi);
      if (id !== undefined && !ringCenter.has(id)) ringCenter.set(id, ri.center);
    }
  }
  return { hCount, badValence, ringCenter };
}

/** Parse SMILES into a document fragment centred on (0,0). Throws on invalid input. */
export function fromSmiles(smiles: string): Doc {
  const prepared = smiles.trim().replace(/\[\*\]|\*/g, '[Ts]');
  const mol = Molecule.fromSmiles(prepared);
  mol.inventCoordinates();
  return fromOCL(mol);
}

export function fromMolfile(text: string): Doc {
  const mol = Molecule.fromMolfile(starToOcl(text));
  if (!mol.getAllAtoms()) throw new Error('No atoms found in MOL file');
  let has2D = false;
  for (let i = 0; i < mol.getAllAtoms(); i++) if (mol.getAtomX(i) !== 0 || mol.getAtomY(i) !== 0) has2D = true;
  if (!has2D) mol.inventCoordinates();
  return fromOCL(mol);
}

export function fromOCL(mol: Mol): Doc {
  const doc = emptyDoc();
  mol.ensureHelperArrays(Molecule.cHelperNeighbours);
  const n = mol.getAllAtoms();
  const scale = BOND / (mol.getAverageBondLength(false) || 1);
  let cx = 0, cy = 0;
  for (let i = 0; i < n; i++) { cx += mol.getAtomX(i); cy += mol.getAtomY(i); }
  cx /= n || 1; cy /= n || 1;
  const ids: number[] = [];
  for (let i = 0; i < n; i++) {
    const z = mol.getAtomicNo(i);
    const el = z === STAR_PLACEHOLDER || z === 0 ? '*' : Molecule.cAtomLabel[z];
    const id = doc.nextId++;
    const atom = { id, el, x: (mol.getAtomX(i) - cx) * scale, y: (mol.getAtomY(i) - cy) * scale, charge: mol.getAtomCharge(i) } as Doc['atoms'][number];
    const mass = mol.getAtomMass(i);
    if (mass) atom.iso = mass;
    doc.atoms.push(atom);
    ids.push(id);
  }
  for (let b = 0; b < mol.getAllBonds(); b++) {
    const type = mol.getBondType(b);
    const order = Math.min(3, Math.max(1, mol.getBondOrder(b))) as 1 | 2 | 3;
    let stereo: 'none' | 'wedge' | 'hash' = 'none';
    if (type === Molecule.cBondTypeUp) stereo = 'wedge';
    else if (type === Molecule.cBondTypeDown) stereo = 'hash';
    doc.bonds.push({ id: doc.nextId++, a: ids[mol.getBondAtom(0, b)], b: ids[mol.getBondAtom(1, b)], order, stereo });
  }
  return doc;
}

export function toSmiles(doc: Doc): string {
  if (!doc.atoms.length) return '';
  const { mol } = toOCL(doc);
  mol.ensureHelperArrays(Molecule.cHelperParities);
  return mol.toIsomericSmiles().replace(/\[\?\]/g, '*');
}

export function toMolfile(doc: Doc, v3000 = false): string {
  const { mol } = toOCL(doc);
  // OCL's internal y axis points down; molfiles expect y up, which toMolfile handles.
  return oclToStar(v3000 ? mol.toMolfileV3() : mol.toMolfile());
}

// OCL writes atoms with atomic number 0 as '?', which most programs reject; '*' is the usual
// symbol for a polymer end / attachment point (read by RDKit, Open Babel, ChemDraw).
const V2000_ATOM = /^(\s*-?\d+\.\d{4}\s+-?\d+\.\d{4}\s+-?\d+\.\d{4}) (\?|\*)( *)(?=\s+-?\d)/gm;
const V3000_ATOM = /^(M {2}V30 \d+) (\?|\*) /gm;
function oclToStar(molfile: string): string {
  return molfile.replace(V2000_ATOM, (_m, xyz, _s, pad) => `${xyz} *${pad}`).replace(V3000_ATOM, '$1 * ');
}
function starToOcl(molfile: string): string {
  return molfile.replace(V2000_ATOM, (_m, xyz, _s, pad) => `${xyz} ?${pad}`).replace(V3000_ATOM, '$1 ? ');
}

/** Relayout with OCL's 2D coordinate generator, keeping the drawing's centre. */
export function cleanCoordinates(doc: Doc): Map<number, Point> {
  const { mol, ids } = toOCL(doc);
  const out = new Map<number, Point>();
  if (!ids.length) return out;
  mol.inventCoordinates();
  const fresh = fromOCL(mol);
  const oldC = centerOf(doc.atoms);
  fresh.atoms.forEach((a, i) => out.set(ids[i], { x: a.x + oldC.x, y: a.y + oldC.y }));
  return out;
}

function centerOf(pts: Point[]): Point {
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const p of pts) { x0 = Math.min(x0, p.x); y0 = Math.min(y0, p.y); x1 = Math.max(x1, p.x); y1 = Math.max(y1, p.y); }
  return { x: (x0 + x1) / 2, y: (y0 + y1) / 2 };
}

// Standard atomic weights (IUPAC 2021, abridged) used for elemental analysis.
const WEIGHTS: Record<string, number> = {
  H: 1.008, He: 4.0026, Li: 6.94, Be: 9.0122, B: 10.81, C: 12.011, N: 14.007, O: 15.999, F: 18.998, Ne: 20.18,
  Na: 22.99, Mg: 24.305, Al: 26.982, Si: 28.085, P: 30.974, S: 32.06, Cl: 35.45, Ar: 39.95, K: 39.098, Ca: 40.078,
  Ti: 47.867, V: 50.942, Cr: 51.996, Mn: 54.938, Fe: 55.845, Co: 58.933, Ni: 58.693, Cu: 63.546, Zn: 65.38,
  Ga: 69.723, Ge: 72.63, As: 74.922, Se: 78.971, Br: 79.904, Rb: 85.468, Sr: 87.62, Zr: 91.224, Mo: 95.95,
  Ru: 101.07, Rh: 102.91, Pd: 106.42, Ag: 107.87, Cd: 112.41, In: 114.82, Sn: 118.71, Sb: 121.76, Te: 127.6,
  I: 126.9, Cs: 132.91, Ba: 137.33, Pt: 195.08, Au: 196.97, Hg: 200.59, Pb: 207.2, Bi: 208.98, Os: 190.23, Ir: 192.22, W: 183.84,
};

export interface Properties {
  formula: string;
  mw: number;
  exactMass: number;
  composition: { el: string; pct: number }[];
  logP?: number;
  tpsa?: number;
  hbd?: number;
  hba?: number;
  rotatable?: number;
  stereocenters?: number;
  isPolymer: boolean;
  fragments: number;
}

export function properties(doc: Doc): Properties | null {
  const real = doc.atoms.filter((a) => a.el !== '*');
  if (!real.length) return null;
  const { mol } = toOCL(doc);
  const mf = mol.getMolecularFormula();
  const counts = parseFormula(mf.formula);
  const composition: { el: string; pct: number }[] = [];
  let total = 0;
  let known = true;
  for (const [el, n] of counts) {
    if (!(el in WEIGHTS)) known = false;
    else total += WEIGHTS[el] * n;
  }
  if (known && total > 0) {
    for (const [el, n] of hillOrder(counts)) composition.push({ el, pct: (100 * WEIGHTS[el] * n) / total });
  }
  const isPolymer = repeatUnitEnds(doc).size > 0;
  const props: Properties = {
    formula: hillFormula(counts, mol),
    mw: mf.relativeWeight,
    exactMass: mf.absoluteWeight,
    composition,
    isPolymer,
    fragments: fragments(doc).length,
  };
  try {
    const p = new MoleculeProperties(mol);
    props.logP = p.logP;
    props.tpsa = p.polarSurfaceArea;
    props.hbd = p.donorCount;
    props.hba = p.acceptorCount;
    props.rotatable = p.rotatableBondCount;
    props.stereocenters = p.stereoCenterCount;
  } catch { /* properties unavailable for unusual structures */ }
  return props;
}

function parseFormula(f: string): Map<string, number> {
  const m = new Map<string, number>();
  const re = /([A-Z][a-z]?)(\d*)/g;
  let r: RegExpExecArray | null;
  while ((r = re.exec(f))) {
    if (!r[1]) continue;
    m.set(r[1], (m.get(r[1]) ?? 0) + (r[2] ? parseInt(r[2], 10) : 1));
  }
  return m;
}

function hillOrder(counts: Map<string, number>): [string, number][] {
  const entries = [...counts.entries()];
  const hasC = counts.has('C');
  return entries.sort((a, b) => {
    const rank = (e: string) => (hasC ? (e === 'C' ? 0 : e === 'H' ? 1 : 2) : 2);
    return rank(a[0]) - rank(b[0]) || a[0].localeCompare(b[0]);
  });
}

function hillFormula(counts: Map<string, number>, mol: Mol): string {
  let charge = 0;
  for (let i = 0; i < mol.getAllAtoms(); i++) charge += mol.getAtomCharge(i);
  const body = hillOrder(counts).map(([el, n]) => (n === 1 ? el : `${el}${n}`)).join('');
  if (!charge) return body;
  const mag = Math.abs(charge);
  return `${body}${mag > 1 ? mag : ''}${charge > 0 ? '+' : '−'}`;
}
