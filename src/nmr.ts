// Predicted ¹H and ¹³C NMR spectra by additive increments.
//
// ¹H, sp³ C–H: base value for CH₃ / CH₂ / CH plus α and β substituent increments
//   (Curphy–Morrison additivity constants, as tabulated in Pavia, Lampman, Kriz & Vyvyan,
//   Introduction to Spectroscopy, 5th ed., Appendix 6), with geminal substituents damped.
// ¹H, alkene: δ = 5.25 + Z_gem + Z_cis + Z_trans (Pascual, Meier & Simon, Helv. Chim. Acta 49 (1966) 164).
// ¹H, benzene: δ = 7.26 + Z_ortho + Z_meta + Z_para; ¹³C benzene 128.5 + Z_ipso/o/m/p;
// ¹³C alkene 123.3 + Z_α + Z_β; substituent tables from E. Pretsch, P. Bühlmann & M. Badertscher,
//   Structure Determination of Organic Compounds, 4th ed., Springer (2009).
// ¹³C, sp³ carbon: Lindeman–Adams equation for the carbon skeleton (Anal. Chem. 43 (1971) 1245) plus
//   α/β/γ substituent increments for the functional groups (Pretsch).
// Couplings are first order (n+1 rule) with typical J values; couplings below 2.5 Hz are ignored.
//
// For a polymer repeat unit (two * ends) the unit is treated as part of an endless chain: the atom next to
// one * continues into the next unit at the other *, so backbone atoms see the neighbours they have in the
// real chain. Polymer signals are drawn broad and without fine structure, as they appear in practice.
//
// These are estimates: typically ±0.2–0.3 ppm for ¹H and ±3–5 ppm for ¹³C. Peaks the tables only cover
// approximately are flagged `rough`; elements and species outside the tables are refused.
import { Molecule } from 'openchemlib';
import { Doc, validBrackets } from './model';
import { toOCL } from './chem';

export type Nucleus = '1H' | '13C';
export const NMR_SOLVENTS = ['CDCl3', 'DMSO-d6', 'CD3OD', 'D2O', 'acetone-d6', 'C6D6', 'CD3CN', 'CD2Cl2'] as const;
export type NmrSolvent = (typeof NMR_SOLVENTS)[number];

export interface Coupling { n: number; j: number }

export interface NmrPeak {
  /** Centre of the signal, ppm. */
  shift: number;
  /** ¹H: number of protons (integral). ¹³C: number of equivalent carbons. */
  count: number;
  /** Atom ids in the drawing that give this signal (the carbon or heteroatom carrying the H for ¹H). */
  atoms: number[];
  assignment: string;
  couplings: Coupling[];
  /** s, d, t, q, quint, sext, sept, dd, dt, td, tt, dq, qd, m, br s, br */
  mult: string;
  /** Line width (FWHM), ppm. Broad for polymers and exchangeable protons. */
  width: number;
  /** Where an exchangeable or solvent-dependent signal usually appears. */
  range?: [number, number];
  exchangeable?: boolean;
  rough?: boolean;
  /** ¹³C: no attached H (weak in a routine spectrum). */
  quaternary?: boolean;
  label: string;
}

export interface NmrSpectrum {
  nucleus: Nucleus;
  peaks: NmrPeak[];
  isPolymer: boolean;
  notes: string[];
}

export type NmrResult = { ok: true; value: NmrSpectrum } | { ok: false; reason: string };

// ---------------------------------------------------------------------------------------------
// Molecular graph. Edges carry an id so paths can avoid stepping straight back; a polymer's two
// * ends become a "virtual" edge between the atoms next to them (offset ±1 = into the next unit).

interface Edge { id: number; to: number; order: number; arom: boolean; off: number; x: number; y: number }

interface Graph {
  mol: InstanceType<typeof Molecule>;
  ids: number[];
  atoms: number[]; // heavy, non-* atoms
  adj: Edge[][];
  nh: number[];
  hPos: Map<number, { x: number; y: number }[]>; // positions of explicitly drawn H atoms
  x: number[];
  y: number[];
  isPolymer: boolean;
  backbone: Set<number>;
  cls: number[]; // symmetry class (graph refinement)
}

const ALLOWED = new Set([1, 6, 7, 8, 9, 14, 16, 17, 35, 53]);

function buildGraph(doc: Doc): Graph | string {
  if (!doc.atoms.some((a) => a.el !== '*')) return 'Nothing drawn.';
  if (validBrackets(doc).length) return 'Bracketed copolymer blocks are not covered yet: draw one repeat unit with two * ends.';
  const { mol, ids } = toOCL(doc);
  mol.ensureHelperArrays(Molecule.cHelperRings);
  const n = mol.getAllAtoms();
  const isH = (i: number) => mol.getAtomicNo(i) === 1 && !mol.getAtomCharge(i) && (mol.getAtomMass(i) === 0 || mol.getAtomMass(i) === 1);
  const isD = (i: number) => mol.getAtomicNo(i) === 1 && mol.getAtomMass(i) === 2;
  const x: number[] = [], y: number[] = [];
  for (let i = 0; i < n; i++) { x.push(mol.getAtomX(i)); y.push(mol.getAtomY(i)); }
  const nh = new Array<number>(n).fill(0);
  const hPos = new Map<number, { x: number; y: number }[]>();
  const adj: Edge[][] = Array.from({ length: n }, () => []);
  for (let i = 0; i < n; i++) {
    const z = mol.getAtomicNo(i);
    if (isH(i) || isD(i) || z === 0) continue;
    if (!ALLOWED.has(z)) return `${Molecule.cAtomLabel[z]} is not in the NMR increment tables.`;
    if (mol.getAtomRadical(i)) return 'Radicals are paramagnetic: their NMR cannot be predicted this way.';
    nh[i] = mol.getImplicitHydrogens(i);
  }
  let eid = 0;
  for (let b = 0; b < mol.getAllBonds(); b++) {
    const a1 = mol.getBondAtom(0, b), a2 = mol.getBondAtom(1, b);
    if ((isH(a1) || isD(a1)) && (isH(a2) || isD(a2))) return 'H₂ has no useful NMR prediction here.';
    if (isD(a1) || isD(a2)) continue; // deuterium: no ¹H signal
    if (isH(a1) || isH(a2)) {
      const [h, a] = isH(a1) ? [a1, a2] : [a2, a1];
      nh[a]++;
      hPos.set(a, [...(hPos.get(a) ?? []), { x: x[h], y: y[h] }]);
      continue;
    }
    const order = mol.getBondOrder(b);
    const arom = mol.isAromaticBond(b);
    adj[a1].push({ id: eid, to: a2, order, arom, off: 0, x: x[a2], y: y[a2] });
    adj[a2].push({ id: eid, to: a1, order, arom, off: 0, x: x[a1], y: y[a1] });
    eid++;
  }
  const heavy: number[] = [];
  for (let i = 0; i < n; i++) if (!isH(i) && !isD(i)) heavy.push(i);
  // One connected structure.
  const seen = new Set<number>([heavy[0]]);
  const stack = [heavy[0]];
  while (stack.length) for (const e of adj[stack.pop()!]) if (!seen.has(e.to)) { seen.add(e.to); stack.push(e.to); }
  if (heavy.some((i) => !seen.has(i))) return 'Draw or select one structure at a time: mixtures and salts are not covered.';

  const stars = heavy.filter((i) => mol.getAtomicNo(i) === 0);
  let isPolymer = false;
  const backbone = new Set<number>();
  if (stars.length) {
    if (stars.length !== 2 || stars.some((s) => adj[s].length !== 1 || adj[s][0].order !== 1)) {
      return 'Open * attachment points: NMR needs a whole molecule or a repeat unit with two * ends.';
    }
    const [s1, s2] = stars;
    const a1 = adj[s1][0].to, a2 = adj[s2][0].to;
    // a1's bond to * continues into the next unit at a2 (and the other way round); the neighbour sits where
    // the * is drawn, which keeps cis/trans geometry right.
    const k1 = adj[a1].findIndex((e) => e.to === s1), k2 = adj[a2].findIndex((e) => e.to === s2);
    const id = eid++;
    adj[a1][k1] = { id, to: a2, order: 1, arom: false, off: -1, x: x[s1], y: y[s1] };
    adj[a2][k2] = { id, to: a1, order: 1, arom: false, off: +1, x: x[s2], y: y[s2] };
    adj[s1] = []; adj[s2] = [];
    isPolymer = true;
    // Backbone: the shortest real path from a1 to a2.
    const prev = new Map<number, number>([[a1, -1]]);
    const q = [a1];
    while (q.length) {
      const i = q.shift()!;
      if (i === a2) break;
      for (const e of adj[i]) if (e.off === 0 && !prev.has(e.to)) { prev.set(e.to, i); q.push(e.to); }
    }
    for (let i: number | undefined = a2; i !== undefined && i >= 0; i = prev.get(i)) backbone.add(i);
  }
  const atoms = heavy.filter((i) => mol.getAtomicNo(i) !== 0);
  for (const i of atoms) {
    const z = mol.getAtomicNo(i), c = mol.getAtomCharge(i);
    if (!c) continue;
    const nitro = z === 7 && c === 1 && adj[i].filter((e) => mol.getAtomicNo(e.to) === 8).length >= 2;
    const nitroO = z === 8 && c === -1 && adj[i].length === 1 && mol.getAtomicNo(adj[i][0].to) === 7 && mol.getAtomCharge(adj[i][0].to) === 1;
    const ammonium = z === 7 && c === 1 && adj[i].every((e) => e.order === 1);
    if (!nitro && !nitroO && !ammonium) return 'Charged species other than nitro groups and ammonium are not covered.';
  }

  // Symmetry classes by iterative refinement over the (multi)graph.
  let cls = heavy.map(() => 0);
  const pos = new Map<number, number>(); heavy.forEach((h, k) => pos.set(h, k));
  const keyOf = (i: number) => `${mol.getAtomicNo(i)}|${nh[i]}|${mol.isAromaticAtom(i) ? 1 : 0}|${mol.getAtomCharge(i)}`;
  const init = heavy.map(keyOf);
  const uniq = [...new Set(init)].sort();
  cls = init.map((k) => uniq.indexOf(k));
  for (let it = 0; it < heavy.length; it++) {
    const sig = heavy.map((h, k) => `${cls[k]}:${adj[h].map((e) => `${e.arom ? 'a' : e.order}-${cls[pos.get(e.to)!]}`).sort().join(',')}`);
    const u = [...new Set(sig)].sort();
    const next = sig.map((s) => u.indexOf(s));
    const done = u.length === new Set(cls).size;
    cls = next;
    if (done) break;
  }
  const clsByAtom = new Array<number>(n).fill(-1);
  heavy.forEach((h, k) => { clsByAtom[h] = cls[k]; });
  return { mol, ids, atoms, adj, nh, hPos, x, y, isPolymer, backbone, cls: clsByAtom };
}

// ---------------------------------------------------------------------------------------------
// Substituent classification: what is attached to atom `c` through edge `e`.

type Kind =
  | 'alkyl' | 'vinyl' | 'alkynyl' | 'aryl' | 'CHO' | 'COR' | 'COOH' | 'COOR' | 'CONR' | 'COCl' | 'CN' | 'imine'
  | 'OH' | 'OR' | 'OAr' | 'OCOR' | 'NH2' | 'NHR' | 'NR2' | 'NAr' | 'NCOR' | 'Nring' | 'Nplus' | 'NO2'
  | 'SH' | 'SR' | 'SO' | 'SO2' | 'F' | 'Cl' | 'Br' | 'I' | 'Si';

class Chem {
  constructor(readonly g: Graph) {}
  z(i: number) { return this.g.mol.getAtomicNo(i); }
  arom(i: number) { return this.g.mol.isAromaticAtom(i); }
  nh(i: number) { return this.g.nh[i]; }
  edges(i: number) { return this.g.adj[i]; }
  /** Edges of `i` other than `e` (an edge object from i's own list). */
  others(i: number, e: Edge) { return this.g.adj[i].filter((f) => f !== e); }
  /** The reverse of edge `e`, i.e. the matching entry in the far atom's list. */
  back(e: Edge, _from?: number) { return this.g.adj[e.to].find((f) => f.id === e.id && f.off === -e.off)!; }
  /** Double bond to a terminal O (carbonyl / sulfonyl O). */
  oxo(i: number) { return this.g.adj[i].filter((e) => e.order === 2 && this.z(e.to) === 8 && this.g.adj[e.to].length === 1); }
  isCarbonyl(i: number) { return this.z(i) === 6 && this.oxo(i).length === 1; }
  /** sp³ carbon: no multiple bonds, not aromatic. */
  sp3C(i: number) { return this.z(i) === 6 && !this.arom(i) && this.g.adj[i].every((e) => e.order === 1); }

  kind(c: number, e: Edge): Kind | null {
    const x = e.to;
    const z = this.z(x);
    const rest = this.others(x, this.back(e, c));
    if (z === 6) {
      if (this.arom(x)) return 'aryl';
      if (rest.some((f) => f.order === 3)) return rest.some((f) => f.order === 3 && this.z(f.to) === 7) ? 'CN' : 'alkynyl';
      if (this.oxo(x).length) {
        const sub = rest.filter((f) => !(f.order === 2 && this.z(f.to) === 8));
        if (sub.some((f) => this.z(f.to) === 8)) return sub.some((f) => this.z(f.to) === 8 && (this.nh(f.to) > 0 || this.g.mol.getAtomCharge(f.to) === -1)) ? 'COOH' : 'COOR';
        if (sub.some((f) => this.z(f.to) === 7)) return 'CONR';
        if (sub.some((f) => this.z(f.to) === 17 || this.z(f.to) === 35)) return 'COCl';
        return this.nh(x) > 0 ? 'CHO' : 'COR';
      }
      if (rest.some((f) => f.order === 2 && this.z(f.to) === 7)) return 'imine';
      if (rest.some((f) => f.order === 2)) return 'vinyl';
      return 'alkyl';
    }
    if (z === 8) {
      if (this.nh(x) > 0) return 'OH';
      if (this.g.mol.getAtomCharge(x) === -1) return null;
      const y = rest[0];
      if (!y) return null;
      if (this.isCarbonyl(y.to) || (this.z(y.to) === 16 && this.oxo(y.to).length)) return 'OCOR';
      if (this.arom(y.to)) return 'OAr';
      return 'OR';
    }
    if (z === 7) {
      if (this.g.mol.getAtomCharge(x) === 1) return this.edges(x).filter((f) => this.z(f.to) === 8).length >= 2 ? 'NO2' : 'Nplus';
      if (this.arom(x)) return 'Nring';
      if (rest.some((f) => f.order > 1)) return 'imine';
      if (rest.some((f) => this.isCarbonyl(f.to) || (this.z(f.to) === 16 && this.oxo(f.to).length))) return 'NCOR';
      if (rest.some((f) => this.arom(f.to))) return 'NAr';
      return this.nh(x) === 2 ? 'NH2' : this.nh(x) === 1 ? 'NHR' : 'NR2';
    }
    if (z === 16) {
      const o = this.oxo(x).length;
      return o >= 2 ? 'SO2' : o === 1 ? 'SO' : this.nh(x) ? 'SH' : 'SR';
    }
    if (z === 9) return 'F';
    if (z === 17) return 'Cl';
    if (z === 35) return 'Br';
    if (z === 53) return 'I';
    if (z === 14) return 'Si';
    return null;
  }

  /** Alkyl carbon subtype for ring and alkene increments. */
  alkylSub(c: number, e: Edge): 'Me' | 'CH2R' | 'CHR2' | 'CR3' | 'CH2X' | 'CF3' {
    const x = e.to;
    const rest = this.others(x, this.back(e, c));
    const nF = rest.filter((f) => this.z(f.to) === 9).length;
    if (nF >= 2) return 'CF3';
    if (rest.some((f) => ![1, 6].includes(this.z(f.to)) || this.arom(f.to) || this.isCarbonyl(f.to))) return 'CH2X';
    const heavyN = rest.length;
    return heavyN === 0 ? 'Me' : heavyN === 1 ? 'CH2R' : heavyN === 2 ? 'CHR2' : 'CR3';
  }
}

// ---------------------------------------------------------------------------------------------
// ¹H tables

// sp³ C–H α increments for [CH₃, CH₂, CH] and β increment (all), Curphy–Morrison style.
const H_SP3: Partial<Record<Kind, [number, number, number, number]>> = {
  alkyl: [0, 0, 0, 0],
  vinyl: [0.78, 0.75, 0.70, 0.10],
  alkynyl: [0.90, 0.90, 0.90, 0.15],
  aryl: [1.45, 1.45, 1.35, 0.35],
  CHO: [1.30, 1.10, 0.95, 0.25],
  COR: [1.20, 1.20, 1.00, 0.15],
  COOH: [1.20, 1.15, 1.00, 0.25],
  COOR: [1.10, 1.10, 0.95, 0.25],
  CONR: [1.10, 1.00, 0.90, 0.20],
  COCl: [1.80, 1.80, 1.60, 0.40],
  CN: [1.10, 1.15, 1.10, 0.40],
  imine: [1.00, 1.00, 1.00, 0.15],
  OH: [2.50, 2.40, 2.30, 0.25],
  OR: [2.40, 2.20, 2.00, 0.25],
  OAr: [2.85, 2.75, 2.85, 0.40],
  OCOR: [2.80, 2.95, 3.40, 0.35],
  NH2: [1.45, 1.40, 1.40, 0.15],
  NHR: [1.45, 1.40, 1.40, 0.15],
  NR2: [1.40, 1.30, 1.30, 0.15],
  NAr: [1.95, 2.00, 2.00, 0.20],
  NCOR: [1.90, 2.05, 2.45, 0.20],
  Nring: [2.70, 2.60, 2.60, 0.35],
  Nplus: [2.30, 2.30, 2.30, 0.40],
  NO2: [3.40, 3.20, 3.10, 0.65],
  SH: [1.20, 1.30, 1.60, 0.30],
  SR: [1.20, 1.30, 1.60, 0.30],
  SO: [1.70, 1.60, 1.60, 0.30],
  SO2: [1.95, 1.85, 1.80, 0.40],
  F: [3.35, 3.15, 3.10, 0.45],
  Cl: [2.15, 2.30, 2.55, 0.55],
  Br: [1.80, 2.20, 2.65, 0.75],
  I: [1.25, 2.00, 2.70, 0.90],
  Si: [-0.90, -0.60, -0.50, 0],
};
const H_BASE = [0.23, 0.90, 1.20, 1.55]; // CH₄ (by H count 4) is handled via index: [CH4, CH3, CH2, CH]

// Alkene H: [gem, cis, trans] (Pascual–Meier–Simon).
const H_ALKENE: Record<string, [number, number, number]> = {
  alkyl: [0.45, -0.22, -0.28],
  CH2X: [0.64, -0.01, -0.02],
  vinyl: [1.00, -0.09, -0.23],
  alkynyl: [0.47, 0.38, 0.12],
  aryl: [1.38, 0.36, -0.07],
  CHO: [1.02, 0.95, 1.17],
  COR: [1.10, 1.12, 0.87],
  COOH: [0.97, 1.41, 0.71],
  COOR: [0.80, 1.18, 0.55],
  CONR: [1.37, 0.98, 0.46],
  COCl: [1.11, 1.46, 1.01],
  CN: [0.27, 0.75, 0.55],
  imine: [0.9, 0.9, 0.6],
  OH: [1.22, -1.07, -1.21],
  OR: [1.22, -1.07, -1.21],
  OAr: [1.21, -0.60, -1.00],
  OCOR: [2.11, -0.35, -0.64],
  NH2: [0.80, -1.26, -1.21],
  NHR: [0.80, -1.26, -1.21],
  NR2: [0.80, -1.26, -1.21],
  NAr: [0.80, -1.26, -1.21],
  NCOR: [2.08, -0.57, -0.72],
  Nring: [2.08, -0.57, -0.72],
  NO2: [1.87, 1.32, 0.62],
  SH: [1.11, -0.29, -0.13],
  SR: [1.11, -0.29, -0.13],
  SO: [1.27, 0.67, 0.41],
  SO2: [1.55, 1.16, 0.93],
  F: [1.54, -0.40, -1.02],
  Cl: [1.08, 0.18, 0.13],
  Br: [1.07, 0.45, 0.55],
  I: [1.14, 0.81, 0.88],
  Si: [0.90, 0.90, 0.60],
};

// Benzene ring H: [ortho, meta, para].
const H_AROM: Record<string, [number, number, number]> = {
  Me: [-0.18, -0.11, -0.21], CH2R: [-0.14, -0.06, -0.17], CHR2: [-0.13, -0.08, -0.18], CR3: [0.02, -0.08, -0.21],
  CH2X: [-0.07, -0.07, -0.07], CF3: [0.32, 0.14, 0.20],
  vinyl: [0.04, -0.04, -0.12], alkynyl: [0.15, -0.02, -0.01], aryl: [0.30, 0.12, 0.10],
  CHO: [0.56, 0.22, 0.29], COR: [0.62, 0.14, 0.21], COOH: [0.85, 0.18, 0.27], COOR: [0.71, 0.11, 0.21],
  CONR: [0.61, 0.10, 0.17], COCl: [0.84, 0.22, 0.36], CN: [0.36, 0.18, 0.28], imine: [0.5, 0.1, 0.1],
  OH: [-0.56, -0.12, -0.45], OR: [-0.48, -0.09, -0.44], OAr: [-0.29, -0.05, -0.23], OCOR: [-0.25, 0.03, -0.13],
  NH2: [-0.75, -0.25, -0.65], NHR: [-0.80, -0.22, -0.68], NR2: [-0.66, -0.18, -0.67], NAr: [-0.55, -0.15, -0.55],
  NCOR: [0.12, -0.07, -0.28], Nring: [0.0, 0.0, 0.0], Nplus: [0.6, 0.3, 0.3], NO2: [0.95, 0.26, 0.38],
  SH: [-0.08, -0.16, -0.22], SR: [-0.08, -0.10, -0.24], SO: [0.43, 0.24, 0.24], SO2: [0.76, 0.35, 0.45],
  F: [-0.26, 0.00, -0.20], Cl: [0.03, -0.02, -0.09], Br: [0.18, -0.08, -0.04], I: [0.39, -0.21, 0.00],
  Si: [0.22, -0.02, -0.02],
};

// ---------------------------------------------------------------------------------------------
// ¹³C tables

// sp³ carbon substituent increments [α, β, γ] (Pretsch, n-alkyl positions).
const C_SP3: Partial<Record<Kind, [number, number, number]>> = {
  vinyl: [21.5, 6.9, -2.1], alkynyl: [4.4, 5.6, -3.4], aryl: [23.0, 9.3, -2.6],
  CHO: [31.0, 0, -2.7], COR: [32.0, 1.0, -2.0], COOH: [22.0, 2.0, -2.8], COOR: [22.5, 2.0, -2.8],
  CONR: [24.0, 2.6, -3.2], COCl: [33.1, 2.3, -3.6], CN: [3.1, 2.4, -3.3], imine: [25, 3, -2],
  OH: [49.0, 10.1, -6.2], OR: [58.0, 7.2, -5.8], OAr: [57.0, 7.0, -5.8], OCOR: [52.0, 6.5, -6.0],
  NH2: [28.3, 11.3, -5.1], NHR: [36.9, 8.1, -4.6], NR2: [40.5, 5.0, -4.6], NAr: [33.0, 7.0, -4.6],
  NCOR: [28.5, 7.0, -5.0], Nring: [38.0, 7.0, -5.0], Nplus: [32.0, 6.0, -5.0], NO2: [61.6, 3.1, -4.6],
  SH: [10.6, 11.4, -3.6], SR: [20.4, 6.2, -2.7], SO: [43.0, 1.0, -3.0], SO2: [46.5, 1.0, -3.0],
  F: [70.1, 7.8, -6.8], Cl: [31.0, 10.0, -5.1], Br: [18.9, 11.0, -3.8], I: [-7.2, 10.9, -1.5],
  Si: [2.3, 1.0, -1.0],
};
// Lindeman–Adams constants.
const LA_B = [-2.3, 6.80, 15.34, 23.46, 27.77];
const LA_A: Record<string, number> = {
  '12': 9.56, '13': 17.83, '14': 25.48, '22': 9.75, '23': 16.70, '24': 21.43,
  '32': 6.60, '33': 11.14, '34': 14.70, '42': 2.26, '43': 3.96, '44': 7.35,
};
const LA_GAMMA = [0, -2.99, -2.69, -2.07, 0.68];
const LA_DELTA = [0, 0.49, 0.48, 0, 0];
const RING_C: Record<number, number> = { 3: -37, 4: -12, 5: -9, 6: -5 };

// Benzene ¹³C: [ipso, ortho, meta, para].
const C_AROM: Record<string, [number, number, number, number]> = {
  Me: [9.3, 0.7, -0.1, -2.9], CH2R: [15.6, -0.5, 0.0, -2.6], CHR2: [20.1, -2.0, 0.0, -2.5], CR3: [22.1, -3.4, -0.4, -3.1],
  CH2X: [12.4, -1.2, 0.2, -1.1], CF3: [2.6, -3.1, 0.4, 3.4],
  vinyl: [8.9, -2.3, -0.1, -0.8], alkynyl: [-6.2, 3.6, -0.4, -0.3], aryl: [13.1, -1.1, 0.4, -1.2],
  CHO: [8.2, 1.2, 0.5, 5.8], COR: [8.9, 0.0, -0.1, 4.4], COOH: [2.1, 1.6, -0.1, 5.2], COOR: [2.0, 1.2, -0.1, 4.3],
  CONR: [5.0, -1.2, 0.1, 3.4], COCl: [4.7, 2.7, 0.3, 6.6], CN: [-16.0, 3.5, 0.7, 4.3], imine: [7.5, -1.0, 0.0, 2.5],
  OH: [26.9, -12.8, 1.4, -7.4], OR: [31.4, -14.4, 1.0, -7.7], OAr: [27.6, -11.2, -0.3, -6.9], OCOR: [22.4, -7.1, 0.4, -3.2],
  NH2: [18.2, -13.4, 0.8, -10.0], NHR: [21.4, -16.2, 0.8, -11.6], NR2: [22.4, -15.7, 0.8, -11.8], NAr: [14.7, -10.6, 0.9, -7.6],
  NCOR: [9.7, -8.1, 0.2, -4.4], Nring: [11.0, -8.0, 1.0, -2.0], Nplus: [19.0, -6.0, 1.5, 1.0], NO2: [19.9, -4.9, 0.9, 6.1],
  SH: [2.1, 0.7, 0.3, -3.2], SR: [10.0, -1.9, 0.2, -3.6], SO: [17.6, -5.0, 1.0, 2.4], SO2: [12.3, -1.1, 0.8, 5.1],
  F: [34.8, -13.0, 1.6, -4.4], Cl: [6.3, 0.4, 1.4, -1.9], Br: [-5.8, 3.2, 1.6, -1.6], I: [-34.1, 8.9, 1.6, -1.1],
  Si: [11.6, 4.9, -0.7, 0.4],
};

// Alkene ¹³C: [same carbon, other carbon]. Alkyl substituents use the Pascual-style α/β/γ terms.
const C_ALKENE: Record<string, [number, number]> = {
  vinyl: [13.6, -7.0], alkynyl: [-6.0, 5.0], aryl: [12.5, -11.0],
  CHO: [13.1, 12.7], COR: [13.8, 4.7], COOH: [4.2, 8.9], COOR: [6.3, 7.0], CONR: [5.0, 5.0], COCl: [8.1, 14.0],
  CN: [-15.1, 14.2], imine: [10, 2],
  OH: [29.4, -38.9], OR: [29.4, -38.9], OAr: [22.0, -32.0], OCOR: [18.4, -26.7],
  NH2: [28.0, -32.0], NHR: [28.0, -32.0], NR2: [28.0, -32.0], NAr: [24.0, -28.0], NCOR: [6.5, -29.2], Nring: [6.5, -29.2],
  Nplus: [19.8, -10.6], NO2: [22.3, -0.9],
  SH: [18.0, -16.0], SR: [18.0, -16.0], SO: [16.0, -1.0], SO2: [14.3, 7.9],
  F: [24.9, -34.3], Cl: [2.6, -6.1], Br: [-8.6, -0.9], I: [-38.1, 7.0], Si: [16.9, 6.7],
};

// ---------------------------------------------------------------------------------------------
// Exchangeable protons: [typical, low, high] by solvent family.

const EXCH: Record<string, { cdcl3: [number, number, number]; dmso: [number, number, number]; label: string }> = {
  OH: { cdcl3: [2.0, 1.0, 5.0], dmso: [4.5, 4.0, 5.5], label: 'O–H (alcohol)' },
  ArOH: { cdcl3: [5.0, 4.0, 7.0], dmso: [9.3, 8.5, 10.0], label: 'O–H (phenol)' },
  COOH: { cdcl3: [11.0, 10.0, 13.0], dmso: [12.2, 11.5, 13.0], label: 'CO–O–H (acid)' },
  NH2: { cdcl3: [1.3, 0.5, 3.0], dmso: [1.5, 0.5, 3.0], label: 'N–H (amine)' },
  ArNH: { cdcl3: [3.6, 3.0, 5.0], dmso: [5.0, 4.0, 6.0], label: 'N–H (aryl amine)' },
  amideNH: { cdcl3: [6.0, 5.0, 8.5], dmso: [8.0, 7.5, 9.0], label: 'N–H (amide)' },
  amideNH2: { cdcl3: [5.8, 5.0, 7.5], dmso: [7.1, 6.8, 7.6], label: 'N–H₂ (amide)' },
  ringNH: { cdcl3: [8.1, 7.5, 11.0], dmso: [11.0, 10.5, 12.0], label: 'N–H (pyrrole, indole)' },
  NplusH: { cdcl3: [7.5, 6.5, 9.0], dmso: [8.0, 7.0, 9.0], label: 'N⁺–H (ammonium)' },
  SH: { cdcl3: [1.4, 1.2, 1.8], dmso: [2.0, 1.5, 2.5], label: 'S–H (thiol)' },
  ArSH: { cdcl3: [3.4, 3.0, 4.0], dmso: [5.0, 4.5, 5.5], label: 'S–H (thiophenol)' },
};

const PERIOD_J = { vic: 7.0, cis: 10.5, trans: 16.5, ortho: 7.8, het5: 3.5, pyr: 5.0, allyl: 6.5 };

// ---------------------------------------------------------------------------------------------

interface HSite {
  atom: number;
  n: number;
  shift: number;
  what: string;
  exch?: { key: string; range: [number, number] };
  rough?: boolean;
  /** Partners: [site index, J Hz, number of H]. Filled after all sites exist. */
  side?: number; // for =CH₂: which side of the C=C the H is on
}

const round = (v: number, d = 2) => Math.round(v * 10 ** d) / 10 ** d;

/** First-order multiplicity label from coupling groups (largest J first). */
function multLabel(groups: Coupling[]): string {
  if (!groups.length) return 's';
  const one = ['s', 'd', 't', 'q', 'quint', 'sext', 'sept'];
  const short = ['', 'd', 't', 'q'];
  if (groups.length === 1) return groups[0].n <= 6 ? one[groups[0].n] : 'm';
  if (groups.length === 2 && groups.every((g) => g.n <= 3)) return groups.map((g) => short[g.n]).join('');
  return 'm';
}

/** Merge couplings with similar J (within 1.2 Hz) into one n+1 group. */
function groupCouplings(list: Coupling[]): Coupling[] {
  const sorted = [...list].filter((c) => c.j >= 2.5 && c.n > 0).sort((a, b) => b.j - a.j);
  const out: Coupling[] = [];
  for (const c of sorted) {
    const g = out.find((o) => Math.abs(o.j - c.j) <= 1.2);
    if (g) { g.j = (g.j * g.n + c.j * c.n) / (g.n + c.n); g.n += c.n; } else out.push({ ...c });
  }
  return out;
}

export function predictNmr(doc: Doc, nucleus: Nucleus, solvent: NmrSolvent = 'CDCl3'): NmrResult {
  let g: Graph | string;
  try { g = buildGraph(doc); } catch { return { ok: false, reason: 'This structure could not be analysed.' }; }
  if (typeof g === 'string') return { ok: false, reason: g };
  return nucleus === '1H' ? predictH(g, solvent) : predictC(g);
}

// ---------------------------------------------------------------------------------------------
// ¹H

function damped(values: number[], weights = [1, 0.75, 0.6, 0.5]): number {
  return [...values].sort((a, b) => Math.abs(b) - Math.abs(a)).reduce((s, v, k) => s + v * (weights[k] ?? 0.5), 0);
}

function ringOf(ch: Chem, c: number, size?: number): number[] | null {
  const rs = ch.g.mol.getRingSet();
  let best: number[] | null = null;
  for (let r = 0; r < rs.getSize(); r++) {
    if (!rs.isAromatic(r)) continue;
    const atoms = rs.getRingAtoms(r);
    if (!atoms.includes(c) || (size && atoms.length !== size)) continue;
    if (!best || atoms.length > best.length) best = atoms; // prefer the 6-ring in fused 5/6 systems
  }
  return best;
}

/** Ring atoms of `ring` ordered starting at `c`, so index = position (0 ipso, 1 ortho …). */
function ringOrder(ch: Chem, ring: number[], c: number): number[] {
  const set = new Set(ring);
  const out = [c];
  let prev = -1, cur = c;
  while (out.length < ring.length) {
    const nx = ch.edges(cur).find((e) => set.has(e.to) && e.to !== prev && !out.includes(e.to));
    if (!nx) break;
    prev = cur; cur = nx.to; out.push(cur);
  }
  return out;
}

interface RingInfo { size: number; hetero: number[]; fused: boolean; order: number[] }

function aromaticRing(ch: Chem, c: number): RingInfo | null {
  const ring = ringOf(ch, c, 6) ?? ringOf(ch, c);
  if (!ring) return null;
  const order = ringOrder(ch, ring, c);
  const hetero = order.map((a, k) => (ch.z(a) !== 6 ? k : -1)).filter((k) => k >= 0);
  const fused = ring.some((a) => isFusion(ch, a, ring));
  return { size: ring.length, hetero, fused, order };
}

/** Ring atom shared with another aromatic ring (naphthalene C4a, indole C3a …). */
function isFusion(ch: Chem, a: number, ring: number[]): boolean {
  return ch.edges(a).some((e) => !ring.includes(e.to) && ch.arom(e.to) && sharesRing(ch, a, e.to));
}

/** Substituents on ring atom `a` (edges leaving the ring, except into a ring fused to this one). */
function ringSubs(ch: Chem, a: number, ring: Set<number>): Edge[] {
  return ch.edges(a).filter((e) => !ring.has(e.to) && !(ch.arom(e.to) && sharesRing(ch, a, e.to)));
}

function sharesRing(ch: Chem, a: number, b: number): boolean {
  const rs = ch.g.mol.getRingSet();
  for (let r = 0; r < rs.getSize(); r++) { const at = rs.getRingAtoms(r); if (at.includes(a) && at.includes(b)) return true; }
  return false;
}

function subKey(ch: Chem, a: number, e: Edge): string | null {
  const k = ch.kind(a, e);
  if (k === 'alkyl') return ch.alkylSub(a, e);
  return k;
}

function predictH(g: Graph, solvent: NmrSolvent): NmrResult {
  const ch = new Chem(g);
  const sites: HSite[] = [];
  const notes: string[] = [];
  const dmsoLike = solvent === 'DMSO-d6';
  const exchanged = solvent === 'D2O' || solvent === 'CD3OD';
  let anyRough = false;

  const exch = (atom: number, n: number, key: string) => {
    const t = EXCH[key];
    const v = dmsoLike ? t.dmso : t.cdcl3;
    sites.push({ atom, n, shift: v[0], what: t.label, exch: { key, range: [v[1], v[2]] } });
  };

  for (const c of g.atoms) {
    const n = g.nh[c];
    if (!n) continue;
    const z = ch.z(c);
    if (z !== 6) {
      const rest = ch.edges(c);
      if (z === 8) {
        if (rest.some((e) => ch.isCarbonyl(e.to))) exch(c, n, 'COOH');
        else if (rest.some((e) => ch.arom(e.to))) exch(c, n, 'ArOH');
        else exch(c, n, 'OH');
      } else if (z === 7) {
        if (g.mol.getAtomCharge(c) === 1) exch(c, n, 'NplusH');
        else if (ch.arom(c)) exch(c, n, 'ringNH');
        else if (rest.some((e) => ch.isCarbonyl(e.to) || (ch.z(e.to) === 16 && ch.oxo(e.to).length))) exch(c, n, n === 2 ? 'amideNH2' : 'amideNH');
        else if (rest.some((e) => ch.arom(e.to))) exch(c, n, 'ArNH');
        else exch(c, n, 'NH2');
      } else if (z === 16) exch(c, n, rest.some((e) => ch.arom(e.to)) ? 'ArSH' : 'SH');
      else if (z === 14) { sites.push({ atom: c, n, shift: 4.2, what: 'Si–H', rough: true }); anyRough = true; }
      continue;
    }
    const edges = ch.edges(c);
    // Aldehyde / formyl H
    if (ch.oxo(c).length) {
      const other = edges.filter((e) => !(e.order === 2 && ch.z(e.to) === 8));
      let s = 9.75, what = 'CHO (aldehyde)';
      if (other.some((e) => ch.z(e.to) === 8)) { s = 8.05; what = 'HCO–O (formate)'; }
      else if (other.some((e) => ch.z(e.to) === 7)) { s = 8.10; what = 'HCO–N (formamide)'; }
      else if (other.some((e) => ch.arom(e.to))) { s = 9.95; what = 'Ar–CHO'; }
      else if (other.some((e) => ch.kind(c, e) === 'vinyl')) { s = 9.50; what = 'C=C–CHO'; }
      sites.push({ atom: c, n, shift: s, what });
      continue;
    }
    // Terminal alkyne
    if (edges.some((e) => e.order === 3)) {
      const t = edges.find((e) => e.order === 3)!;
      const far = ch.others(t.to, ch.back(t, c))[0];
      const k = far ? ch.kind(t.to, far) : null;
      const s = k === 'aryl' ? 3.05 : k && k.startsWith('CO') ? 2.90 : 1.95;
      sites.push({ atom: c, n, shift: s, what: '≡C–H' });
      continue;
    }
    if (ch.arom(c)) {
      const r = aromaticRing(ch, c);
      if (!r) { sites.push({ atom: c, n, shift: 7.2, what: 'Ar–H', rough: true }); anyRough = true; continue; }
      let s = 7.26, rough = r.fused;
      if (r.size === 6) {
        for (const k of r.hetero) {
          const p = Math.min(k, 6 - k);
          s += p === 1 ? 1.34 : p === 3 ? 0.38 : 0;
        }
        if (r.hetero.length > 1) rough = true;
        if (r.fused) {
          // Naphthalene-like: α to the ring fusion 7.84, β 7.47.
          const nextToFusion = [r.order[1], r.order[5]].some((a) => isFusion(ch, a, r.order));
          s += nextToFusion ? 0.58 : 0.21;
        }
      } else if (r.size === 5) {
        rough = rough || r.hetero.length !== 1;
        const h = r.hetero[0] ?? 0;
        const het = ch.z(r.order[h]);
        const alpha = Math.min(h, 5 - h) === 1;
        const base: Record<number, [number, number]> = { 8: [7.42, 6.38], 16: [7.31, 7.10], 7: [6.74, 6.24] };
        const bb = base[het] ?? [7.2, 7.0];
        s = alpha ? bb[0] : bb[1];
      } else { rough = true; s = 7.2; }
      const ring = new Set(r.order);
      for (let k = 1; k < r.order.length; k++) {
        const p = Math.min(k, r.order.length - k); // 1 ortho, 2 meta, 3 para
        for (const e of ringSubs(ch, r.order[k], ring)) {
          const key = subKey(ch, r.order[k], e);
          const z3 = key ? H_AROM[key] : undefined;
          if (!z3) { rough = true; continue; }
          s += z3[Math.min(p, 3) - 1];
        }
      }
      if (r.size === 5 && r.order.some((a, k) => k > 0 && ringSubs(ch, a, ring).length)) rough = true;
      // Pendant rings on a vinyl polymer backbone stack against each other: ortho H ≈ 6.6, meta/para ≈ 7.05 in polystyrene.
      for (let k = 1; k < r.order.length; k++) {
        const p = Math.min(k, r.order.length - k);
        if (ringSubs(ch, r.order[k], ring).some((e) => g.backbone.has(e.to) && pendantAryl(g, ch, e.to))) s -= p === 1 ? 0.5 : 0.1;
      }
      anyRough ||= rough;
      sites.push({ atom: c, n, shift: s, what: r.hetero.length ? 'Het–H (aromatic)' : 'Ar–H', rough });
      continue;
    }
    const dbl = edges.find((e) => e.order === 2);
    if (dbl) {
      if (ch.z(dbl.to) !== 6 || edges.filter((e) => e.order === 2).length > 1 || ch.arom(dbl.to)) {
        sites.push({ atom: c, n, shift: ch.z(dbl.to) === 7 ? 8.1 : 4.7, what: ch.z(dbl.to) === 7 ? 'H–C=N' : '=C–H', rough: true });
        anyRough = true;
        continue;
      }
      const d = dbl.to;
      // Sides are taken against the C=C axis drawn from the lower to the higher atom index, so both carbons of
      // the double bond use the same reference and their H atoms can be compared (cis vs trans coupling).
      const [p0, p1] = c < d ? [c, d] : [d, c];
      const vx = g.x[p1] - g.x[p0], vy = g.y[p1] - g.y[p0];
      const side = (px: number, py: number) => { const v = vx * (py - g.y[p0]) - vy * (px - g.x[p0]); return Math.abs(v) < 1e-3 ? 0 : Math.sign(v); };
      const subsC = edges.filter((e) => e !== dbl);
      const subsD = ch.others(d, ch.back(dbl, c));
      // Sides of this carbon's H atoms.
      let hSides: number[];
      const drawn = g.hPos.get(c);
      if (drawn && drawn.length === n) hSides = drawn.map((p) => side(p.x, p.y));
      else if (n === 2) hSides = [1, -1];
      else hSides = [subsC.length ? -side(subsC[0].x, subsC[0].y) : 0];
      let rough = false;
      const zOf = (atom: number, e: Edge) => {
        const k = subKey(ch, atom, e);
        const key = k === 'Me' || k === 'CH2R' || k === 'CHR2' || k === 'CR3' ? 'alkyl' : k;
        const t = key ? H_ALKENE[key] : undefined;
        if (!t) rough = true;
        return t ?? [0, 0, 0];
      };
      let gem = 0;
      for (const e of subsC) gem += zOf(c, e)[0];
      const dSubs = subsD.map((e) => ({ t: zOf(d, e), s: side(e.x, e.y) }));
      hSides.forEach((hs) => {
        let s = 5.25 + gem;
        for (const ds of dSubs) {
          if (hs === 0 || ds.s === 0) s += (ds.t[1] + ds.t[2]) / 2;
          else s += hs === ds.s ? ds.t[1] : ds.t[2];
        }
        anyRough ||= rough;
        sites.push({ atom: c, n: 1, shift: s, what: '=C–H (alkene)', rough, side: hs });
      });
      continue;
    }
    // sp³ C–H
    const kinds = edges.map((e) => ({ e, k: ch.kind(c, e) }));
    let rough = false;
    const alphas: number[] = [];
    const idx = n >= 4 ? 0 : 4 - n; // CH4 → 0, CH3 → 1, CH2 → 2, CH → 3
    let s = H_BASE[idx];
    for (const { e, k } of kinds) {
      if (!k) { rough = true; continue; }
      const t = H_SP3[k];
      if (!t) { rough = true; continue; }
      if (k !== 'alkyl') {
        let a = t[Math.min(Math.max(idx, 1), 3) - 1];
        // Aryl ketones and aryl esters deshield more than their aliphatic analogues (acetophenone CH₃ 2.6, methyl benzoate OCH₃ 3.9).
        if (k === 'COR' && ch.edges(e.to).some((f) => ch.arom(f.to))) a += 0.45;
        if (k === 'OCOR' && ch.others(e.to, ch.back(e)).some((f) => ch.edges(f.to).some((h) => ch.arom(h.to)))) a += 0.2;
        // Vinyl polymers: the neighbouring pendant rings shield the backbone (polystyrene CH 1.85, not ~2.9).
        if (k === 'aryl' && pendantAryl(g, ch, c)) a -= 1.0;
        alphas.push(a);
      }
      else {
        // β substituents through the alkyl neighbour, damped per β carbon.
        const betas: number[] = [];
        for (const f of ch.others(e.to, ch.back(e, c))) {
          const kb = ch.kind(e.to, f);
          if (!kb || kb === 'alkyl') continue;
          const tb = H_SP3[kb];
          if (!tb) { rough = true; continue; }
          betas.push(kb === 'aryl' && pendantAryl(g, ch, e.to) ? tb[3] * 0.35 : tb[3]);
        }
        s += damped(betas, [1, 0.3, 0.15]);
      }
    }
    if (n >= 4) s = 0.23;
    s += damped(alphas);
    // Geminal fluorines: CHF₂ and CH₂F are outside simple additivity.
    if (edges.filter((e) => ch.z(e.to) === 9).length >= 2) rough = true;
    const ringSize = g.mol.getAtomRingSize(c);
    if (ringSize === 3) s -= 1.0; else if (ringSize === 4) s += 0.75; else if (ringSize === 5 || ringSize === 6) s += alphas.length ? 0 : 0.25;
    anyRough ||= rough;
    const label = n === 3 ? 'CH₃' : n === 2 ? 'CH₂' : n === 1 ? 'CH' : 'CH₄';
    sites.push({ atom: c, n, shift: s, what: describeSp3(g, c, label, kinds.map((k) => k.k)), rough });
  }

  // Couplings between sites on neighbouring atoms (3J), with equivalent sites grouped.
  const groupKey = (s: HSite) => `${g.cls[s.atom]}|${round(s.shift, 2)}|${s.exch ? 'x' : ''}`;
  const groups = new Map<string, number[]>();
  sites.forEach((s, k) => { const key = groupKey(s); groups.set(key, [...(groups.get(key) ?? []), k]); });
  const groupOf = new Map<number, string>();
  for (const [key, list] of groups) for (const k of list) groupOf.set(k, key);
  const sitesOn = new Map<number, number[]>();
  sites.forEach((s, k) => sitesOn.set(s.atom, [...(sitesOn.get(s.atom) ?? []), k]));

  const couplingsOf = (k: number): Coupling[] => {
    const s = sites[k];
    if (s.exch) return [];
    const out: Coupling[] = [];
    const mine = groupOf.get(k);
    const aromC = ch.arom(s.atom);
    for (const e of ch.edges(s.atom)) {
      for (const p of sitesOn.get(e.to) ?? []) {
        const t = sites[p];
        if (t.exch || groupOf.get(p) === mine) continue;
        let j = 0;
        const tArom = ch.arom(t.atom);
        if (aromC && tArom) {
          const r = aromaticRing(ch, s.atom);
          j = r?.size === 5 ? PERIOD_J.het5 : r?.hetero.some((h) => h === 1 || h === r.size - 1) ? PERIOD_J.pyr : PERIOD_J.ortho;
        } else if (aromC || tArom) continue;
        else if (e.order === 2) {
          if (s.side === undefined || t.side === undefined || s.side === 0 || t.side === 0) j = 13;
          else j = s.side === t.side ? PERIOD_J.cis : PERIOD_J.trans;
        } else if (ch.oxo(s.atom).length || ch.oxo(t.atom).length) j = 2.0; // aldehyde H to α-CH
        else if (s.side !== undefined && t.side !== undefined) j = 10.5; // =CH–CH= in a diene
        else if (s.side !== undefined || t.side !== undefined) j = PERIOD_J.allyl;
        else j = PERIOD_J.vic;
        out.push({ n: t.n, j });
      }
    }
    return groupCouplings(out);
  };

  const peaks: NmrPeak[] = [];
  for (const [, list] of groups) {
    const first = sites[list[0]];
    if (first.exch && exchanged) continue;
    const count = list.reduce((s, k) => s + sites[k].n, 0);
    const atoms = [...new Set(list.map((k) => g.ids[sites[k].atom]))];
    const onBackbone = list.some((k) => g.backbone.has(sites[k].atom));
    const couplings = g.isPolymer ? [] : couplingsOf(list[0]);
    let mult = first.exch ? 'br s' : multLabel(couplings);
    let width = first.exch ? 0.08 : 0.0025;
    if (g.isPolymer && !first.exch) { mult = 'br'; width = onBackbone || ch.arom(first.atom) ? 0.18 : 0.06; }
    peaks.push({
      shift: round(first.shift, 2), count, atoms, assignment: first.what, couplings, mult, width,
      range: first.exch?.range, exchangeable: !!first.exch, rough: list.some((k) => sites[k].rough), label: '',
    });
  }
  peaks.sort((a, b) => b.shift - a.shift);
  peaks.forEach((p, k) => { p.label = letter(k); });
  if (exchanged && sites.some((s) => s.exch)) notes.push(`O–H / N–H protons exchange with ${solvent} and are not shown.`);
  if (sites.some((s) => s.exch) && !exchanged) notes.push('O–H, N–H and S–H positions depend on concentration, temperature and water; they are drawn broad.');
  if (anyRough) notes.push('Peaks marked ~ are outside the main increment tables and less reliable.');
  return { ok: true, value: { nucleus: '1H', peaks, isPolymer: g.isPolymer, notes } };
}

/** Backbone sp³ carbon of a vinyl-type polymer chain (both chain neighbours sp³ carbons) carrying an aromatic ring. */
function pendantAryl(g: Graph, ch: Chem, c: number): boolean {
  if (!g.isPolymer || !g.backbone.has(c) || !ch.sp3C(c)) return false;
  const chain = ch.edges(c).filter((e) => g.backbone.has(e.to) || e.off !== 0);
  return chain.length === 2 && chain.every((e) => ch.sp3C(e.to)) && ch.edges(c).some((e) => ch.arom(e.to) && !g.backbone.has(e.to));
}

function letter(k: number): string {
  let s = '';
  k++;
  while (k > 0) { k--; s = String.fromCharCode(97 + (k % 26)) + s; k = Math.floor(k / 26); }
  return s;
}

const KIND_LABEL: Partial<Record<Kind, string>> = {
  aryl: 'Ar', vinyl: 'C=C', alkynyl: 'C≡C', CHO: 'CHO', COR: 'C(=O)R', COOH: 'COOH', COOR: 'C(=O)O', CONR: 'C(=O)N',
  COCl: 'C(=O)Cl', CN: 'C≡N', imine: 'C=N', OH: 'OH', OR: 'O', OAr: 'OAr', OCOR: 'OC(=O)', NH2: 'NH₂', NHR: 'NH', NR2: 'N',
  NAr: 'NAr', NCOR: 'NC(=O)', Nring: 'N(ring)', Nplus: 'N⁺', NO2: 'NO₂', SH: 'SH', SR: 'S', SO: 'S(=O)', SO2: 'SO₂',
  F: 'F', Cl: 'Cl', Br: 'Br', I: 'I', Si: 'Si',
};
// How the group reads next to the observed CHₙ: "Ar–CH₂", "CH₃–OC(=O)" (methyl ester), "CH₂–Cl".
const BEFORE = new Set<Kind>(['aryl', 'vinyl', 'alkynyl']);

function describeSp3(g: Graph, c: number, label: string, kinds: (Kind | null)[]): string {
  const subs = [...new Set(kinds.filter((k): k is Kind => !!k && k !== 'alkyl'))]
    .sort((a, b) => (H_SP3[b]?.[1] ?? 0) - (H_SP3[a]?.[1] ?? 0));
  let text = label;
  if (subs.length) {
    const k = subs[0], t = KIND_LABEL[k] ?? k;
    text = BEFORE.has(k) ? `${t}–${label}` : `${label}–${t}`;
    if (subs.length > 1) text += ` (${subs.slice(1).map((x) => KIND_LABEL[x] ?? x).join(', ')})`;
  }
  return g.isPolymer && g.backbone.has(c) ? `${text}, backbone` : text;
}

// ---------------------------------------------------------------------------------------------
// ¹³C

/** A position on a path: atom index plus which repeat unit it is in (polymers), so paths never revisit an atom. */
interface Step { atom: number; off: number }

function predictC(g: Graph): NmrResult {
  const ch = new Chem(g);
  const notes: string[] = [];
  let anyRough = false;
  const sites: { atom: number; shift: number; what: string; rough: boolean; quaternary: boolean }[] = [];

  const skel = (i: number) => ch.sp3C(i);
  const nSkel = (i: number) => ch.edges(i).filter((e) => skel(e.to)).length;

  for (const c of g.atoms) {
    if (ch.z(c) !== 6) continue;
    const edges = ch.edges(c);
    const nh = g.nh[c];
    let s = 0, what = '', rough = false;
    const nF = edges.filter((e) => ch.z(e.to) === 9).length;
    if (ch.oxo(c).length && !ch.arom(c)) {
      ({ s, what, rough } = carbonyl(ch, c));
    } else if (edges.some((e) => e.order === 3)) {
      const t = edges.find((e) => e.order === 3)!;
      if (ch.z(t.to) === 7) {
        const other = edges.find((e) => e !== t);
        const k = other ? ch.kind(c, other) : null;
        s = k === 'aryl' ? 118.8 : k === 'vinyl' ? 117.5 : 118.5; what = 'C≡N';
      } else {
        const other = edges.find((e) => e !== t);
        const k = other ? ch.kind(c, other) : null;
        if (nh) {
          const far = ch.others(t.to, ch.back(t, c))[0];
          s = far && ch.kind(t.to, far) === 'aryl' ? 77.2 : 68.5;
        } else s = k === 'aryl' ? 89.0 : k === 'alkyl' ? 82.0 : 80.0;
        if (k && k !== 'aryl' && k !== 'alkyl') rough = true;
        what = '–C≡';
      }
    } else if (edges.some((e) => e.order === 2 && ch.z(e.to) === 7) && !ch.arom(c)) {
      s = 160; what = 'C=N'; rough = true;
    } else if (ch.arom(c)) {
      ({ s, what, rough } = aromaticC(ch, c));
    } else if (edges.some((e) => e.order === 2)) {
      ({ s, what, rough } = alkeneC(ch, c));
    } else if (nF >= 2) {
      s = nF >= 3 ? 124 : nh ? 115 : 117; what = nF >= 3 ? 'CF₃' : nh ? 'CHF₂' : 'CF₂'; rough = true;
    } else {
      // Lindeman–Adams skeleton term.
      const sN = Math.min(nSkel(c), 4);
      s = LA_B[sN];
      for (const e of edges) if (skel(e.to)) { const m = Math.min(nSkel(e.to), 4); if (m >= 2) s += LA_A[`${sN}${m}`] ?? 0; }
      // Skeleton carbons by shortest distance (γ = 3 bonds, δ = 4), and β / γ substituents hanging off the
      // skeleton carbons at distance 1 and 2. States carry the repeat-unit offset so a polymer chain is endless.
      const key = (st: Step) => `${st.atom}:${st.off}`;
      const dist = new Map<string, number>([[key({ atom: c, off: 0 }), 0]]);
      let layer: Step[] = [{ atom: c, off: 0 }];
      for (let d = 1; d <= 4 && layer.length; d++) {
        const next: Step[] = [];
        for (const cur of layer) {
          for (const e of ch.edges(cur.atom)) {
            const nx = { atom: e.to, off: cur.off + e.off };
            if (skel(e.to)) {
              if (dist.has(key(nx))) continue;
              dist.set(key(nx), d);
              next.push(nx);
              if (d === 3) s += sN ? LA_GAMMA[sN] : -2.5;
              if (d === 4) s += LA_DELTA[sN] ?? 0;
            } else {
              if (dist.has(key(nx))) continue; // a ring heteroatom reached again the long way round
              dist.set(key(nx), d);
              if (d !== 2 && d !== 3) continue;
              const k = ch.kind(cur.atom, e);
              if (!k || k === 'alkyl') continue;
              const t = C_SP3[k];
              if (!t) { rough = true; continue; }
              s += d === 2 ? t[1] : t[2];
            }
          }
        }
        layer = next;
      }
      // α substituents on this carbon, damped (fluorine strongly).
      const alphas: { v: number; f: boolean; cgroup: boolean }[] = [];
      for (const e of edges) {
        if (skel(e.to)) continue;
        const k = ch.kind(c, e);
        const t = k ? C_SP3[k] : undefined;
        if (!t) { rough = true; continue; }
        // N of a tertiary amide (DMF, PVP): N–C about 7 ppm further downfield than in a secondary amide.
        const v = k === 'NCOR' && ch.nh(e.to) === 0 ? t[0] + 7 : t[0];
        alphas.push({ v, f: k === 'F', cgroup: ch.z(e.to) === 6 });
      }
      const wF = [1, 0.5, 0.12, 0.03], wX = [1, 0.8, 0.6, 0.5];
      // Carbon-type groups add less at branched carbons (Pretsch's n / iso / tert increments: phenyl 23 → 19 → 11).
      const branch = edges.length >= 4 ? 0.65 : edges.length === 3 ? 0.85 : 1;
      for (const a of alphas) if (a.cgroup) a.v *= branch;
      alphas.sort((a, b) => b.v - a.v).forEach((a, k) => { s += a.v * (a.f ? wF[k] : wX[k]) ; });
      const rsz = g.mol.getAtomRingSize(c);
      if (RING_C[rsz]) s += RING_C[rsz];
      const ks = edges.map((e) => ch.kind(c, e));
      what = describeSp3(g, c, nh === 3 ? 'CH₃' : nh === 2 ? 'CH₂' : nh === 1 ? 'CH' : nh === 4 ? 'CH₄' : 'C', ks);
      if (edges.some((e) => ch.edges(e.to).filter((f) => ch.z(f.to) === 9).length >= 2)) rough = true;
    }
    anyRough ||= rough;
    sites.push({ atom: c, shift: s, what, rough, quaternary: nh === 0 });
  }

  const groups = new Map<string, number[]>();
  sites.forEach((s, k) => { const key = `${g.cls[s.atom]}|${round(s.shift, 1)}`; groups.set(key, [...(groups.get(key) ?? []), k]); });
  const peaks: NmrPeak[] = [];
  for (const [, list] of groups) {
    const f = sites[list[0]];
    const backbone = list.some((k) => g.backbone.has(sites[k].atom));
    peaks.push({
      shift: round(Math.max(f.shift, -5), 1), count: list.length, atoms: list.map((k) => g.ids[sites[k].atom]),
      assignment: f.what, couplings: [], mult: g.isPolymer ? 'br' : 's',
      width: g.isPolymer ? (backbone ? 1.2 : 0.5) : 0.02, rough: list.some((k) => sites[k].rough), quaternary: f.quaternary, label: '',
    });
  }
  peaks.sort((a, b) => b.shift - a.shift);
  peaks.forEach((p, k) => { p.label = String(k + 1); });
  if (anyRough) notes.push('Peaks marked ~ are outside the main increment tables and less reliable.');
  return { ok: true, value: { nucleus: '13C', peaks, isPolymer: g.isPolymer, notes } };
}

/** Number of carbons on the α carbon(s) next to a carbonyl, for the branching shift of C=O. */
function alphaBranch(ch: Chem, c: number, e: Edge): number {
  if (!ch.sp3C(e.to)) return 0;
  return ch.others(e.to, ch.back(e, c)).filter((f) => ch.z(f.to) === 6).length;
}

function carbonyl(ch: Chem, c: number): { s: number; what: string; rough: boolean } {
  const subs = ch.edges(c).filter((e) => !(e.order === 2 && ch.z(e.to) === 8));
  const O = subs.filter((e) => ch.z(e.to) === 8), N = subs.filter((e) => ch.z(e.to) === 7);
  const C = subs.filter((e) => ch.z(e.to) === 6), X = subs.filter((e) => [17, 35].includes(ch.z(e.to)));
  const hasH = ch.nh(c) > 0;
  const cType = (e: Edge) => (ch.arom(e.to) ? 'aryl' : ch.kind(c, e) === 'vinyl' ? 'vinyl' : 'alkyl');
  const br = (e: Edge) => [0, 3.5, 6.3, 7.9][Math.min(alphaBranch(ch, c, e), 3)];
  if (O.length === 2) return { s: 155.0, what: 'O–C(=O)–O (carbonate)', rough: false };
  if (O.length === 1 && N.length === 1) return { s: 156.5, what: 'N–C(=O)–O (carbamate)', rough: false };
  if (N.length === 2) return { s: 158.0, what: 'N–C(=O)–N (urea)', rough: false };
  if (O.length === 1) {
    const o = O[0];
    const acid = ch.nh(o.to) > 0 || ch.g.mol.getAtomCharge(o.to) === -1;
    const anhydride = !acid && ch.others(o.to, ch.back(o, c)).some((f) => ch.isCarbonyl(f.to));
    if (hasH) return { s: acid ? 166.0 : 161.0, what: acid ? 'HCOOH' : 'HC(=O)O (formate)', rough: false };
    const t = C[0] ? cType(C[0]) : 'alkyl';
    if (anhydride) return { s: t === 'alkyl' ? 167.0 : 162.5, what: 'C(=O)–O–C(=O) (anhydride)', rough: true };
    if (acid) return { s: t === 'alkyl' ? 177.0 + (C[0] ? br(C[0]) : 0) : 172.0, what: 'COOH', rough: false };
    return { s: t === 'alkyl' ? 171.0 + (C[0] ? br(C[0]) : 0) : t === 'aryl' ? 166.5 : 167.0, what: 'C(=O)O (ester)', rough: false };
  }
  if (N.length === 1) {
    if (hasH) return { s: 162.5, what: 'HC(=O)N (formamide)', rough: false };
    const t = C[0] ? cType(C[0]) : 'alkyl';
    return { s: t === 'alkyl' ? 172.0 + (C[0] ? br(C[0]) / 2 : 0) : t === 'aryl' ? 168.0 : 166.0, what: 'C(=O)N (amide)', rough: false };
  }
  if (X.length) return { s: C[0] && cType(C[0]) === 'aryl' ? 168.0 : 170.0, what: 'C(=O)Cl (acid halide)', rough: false };
  if (hasH) {
    const t = C[0] ? cType(C[0]) : 'alkyl';
    return { s: t === 'aryl' ? 192.0 : t === 'vinyl' ? 193.5 : 200.5 + (C[0] ? br(C[0]) : 0), what: 'CHO', rough: !C[0] };
  }
  if (C.length === 2) {
    const ts = C.map(cType);
    const ring = ch.g.mol.getAtomRingSize(c);
    if (ts.every((t) => t === 'alkyl')) {
      const b = C.reduce((s, e) => s + Math.min(alphaBranch(ch, c, e), 3), 0);
      return { s: 206.7 + 2.3 * b + (ring === 5 ? 6 : ring === 4 ? -2 : 0), what: 'C(=O) (ketone)', rough: ring === 4 || ring === 3 };
    }
    if (ts.filter((t) => t === 'aryl').length === 2) return { s: 196.5, what: 'Ar–C(=O)–Ar', rough: false };
    if (ts.includes('aryl') && ts.includes('alkyl')) return { s: 198.0, what: 'Ar–C(=O)–R', rough: false };
    return { s: 199.0, what: 'C=C–C(=O) (enone)', rough: false };
  }
  return { s: 185, what: 'C=O', rough: true };
}

function aromaticC(ch: Chem, c: number): { s: number; what: string; rough: boolean } {
  const r = aromaticRing(ch, c);
  if (!r) return { s: 128, what: 'Ar C', rough: true };
  let s = 128.5, rough = r.fused;
  if (r.size === 6) {
    for (const k of r.hetero) { const p = Math.min(k, 6 - k); s += k === 0 ? 0 : p === 1 ? 21.4 : p === 2 ? -4.7 : 7.4; }
    if (r.hetero.length > 1) rough = true;
    if (r.fused) {
      const fusion = isFusion(ch, c, r.order);
      const nextToFusion = [r.order[1], r.order[5]].some((a) => isFusion(ch, a, r.order));
      s = fusion ? 133.6 : nextToFusion ? 128.0 : 125.9;
    }
  } else if (r.size === 5) {
    rough = true;
    const h = r.hetero[0] ?? 0;
    const het = ch.z(r.order[h]);
    const alpha = Math.min(h, 5 - h) === 1;
    const base: Record<number, [number, number]> = { 8: [142.8, 109.7], 16: [125.6, 127.3], 7: [118.2, 108.2] };
    s = (base[het] ?? [125, 125])[alpha ? 0 : 1];
  } else rough = true;
  const ring = new Set(r.order);
  let what = 'Ar C–H';
  for (let k = 0; k < r.order.length; k++) {
    const p = Math.min(k, r.order.length - k);
    for (const e of ringSubs(ch, r.order[k], ring)) {
      const key = subKey(ch, r.order[k], e);
      const t = key ? C_AROM[key] : undefined;
      if (!t) { rough = true; continue; }
      if (r.size === 6) s += t[Math.min(p, 3)];
      else if (k === 0) s += t[0] * 0.8;
      if (k === 0) what = `Ar C–${key && ['Me', 'CH2R', 'CHR2', 'CR3', 'CH2X', 'CF3'].includes(key) ? 'C' : KIND_LABEL[key as Kind] ?? key} (ipso)`;
    }
  }
  if (ch.nh(c) === 0 && what === 'Ar C–H') what = 'Ar C (fusion)';
  return { s, what, rough };
}

function alkeneC(ch: Chem, c: number): { s: number; what: string; rough: boolean } {
  const edges = ch.edges(c);
  const dbls = edges.filter((e) => e.order === 2);
  if (dbls.length > 1 || ch.z(dbls[0].to) !== 6) return { s: dbls.length > 1 ? 210 : 140, what: dbls.length > 1 ? '=C= (allene)' : 'C=X', rough: true };
  const dbl = dbls[0];
  const d = dbl.to;
  let s = 123.3, rough = false;
  const alkylTerms = (atom: number, e: Edge, same: boolean) => {
    // α carbon, then β and γ skeleton carbons along the chain.
    s += same ? 10.6 : -7.9;
    const b = ch.others(e.to, ch.back(e, atom)).filter((f) => ch.sp3C(f.to));
    s += b.length * (same ? 7.2 : -1.8);
    for (const f of b) s += ch.others(f.to, ch.back(f, e.to)).filter((h) => ch.sp3C(h.to)).length * (same ? -1.5 : 1.5);
  };
  const apply = (atom: number, list: Edge[], same: boolean) => {
    for (const e of list) {
      const k = ch.kind(atom, e);
      if (k === 'alkyl') { alkylTerms(atom, e, same); continue; }
      const t = k ? C_ALKENE[k] : undefined;
      if (!t) { rough = true; continue; }
      s += same ? t[0] : t[1];
    }
  };
  apply(c, edges.filter((e) => e !== dbl), true);
  apply(d, ch.others(d, ch.back(dbl, c)), false);
  const nh = ch.nh(c);
  return { s, what: nh === 2 ? '=CH₂' : nh === 1 ? '=CH–' : '=C<', rough };
}

// ---------------------------------------------------------------------------------------------
// Spectrum simulation

export interface Line { ppm: number; height: number; width: number; peak: number }

/** Expand peaks into Lorentzian lines (first-order multiplets at `mhz` for ¹H). Heights are per unit area. */
export function spectrumLines(sp: NmrSpectrum, mhz: number): Line[] {
  const lines: Line[] = [];
  const f = sp.nucleus === '1H' ? mhz : mhz * 0.2515;
  sp.peaks.forEach((p, k) => {
    const area = sp.nucleus === '1H' ? p.count : p.count * (p.quaternary ? 0.45 : 1);
    let comps: { off: number; w: number }[] = [{ off: 0, w: 1 }];
    for (const cpl of p.couplings) {
      const next: { off: number; w: number }[] = [];
      const bin = binomial(cpl.n);
      for (const c of comps) bin.forEach((w, i) => next.push({ off: c.off + (i - cpl.n / 2) * (cpl.j / f), w: c.w * w }));
      comps = next;
    }
    const tot = comps.reduce((s, c) => s + c.w, 0);
    for (const c of comps) lines.push({ ppm: p.shift + c.off, height: (area * c.w) / tot, width: p.width, peak: k });
  });
  return lines;
}

function binomial(n: number): number[] {
  let row = [1];
  for (let i = 0; i < n; i++) row = [...row, 0].map((v, k) => v + (k ? row[k - 1] : 0));
  return row;
}

/** Sum of Lorentzians sampled at `n` points from x0 down to x1 (ppm). Lines narrower than `minWidth` are widened for display. */
export function sampleSpectrum(lines: Line[], x0: number, x1: number, n: number, minWidth: number): number[] {
  const ys = new Array<number>(n).fill(0);
  for (const l of lines) {
    const w = Math.max(l.width, minWidth) / 2;
    const amp = l.height / w; // area-normalised: height ∝ area / width
    const reach = w * 60;
    for (let i = 0; i < n; i++) {
      const x = x0 - ((x0 - x1) * i) / (n - 1);
      const dx = x - l.ppm;
      if (Math.abs(dx) > reach) continue;
      ys[i] += amp * (w * w) / (dx * dx + w * w);
    }
  }
  return ys;
}

export function formatShift(p: NmrPeak, nucleus: Nucleus): string {
  const d = nucleus === '1H' ? 2 : 1;
  if (p.range) return `${p.range[0].toFixed(1)}–${p.range[1].toFixed(1)}`;
  if (nucleus === '1H' && p.mult === 'br' && p.width >= 0.1) return `${(p.shift - p.width / 2).toFixed(d - 1)}–${(p.shift + p.width / 2).toFixed(d - 1)}`;
  return p.shift.toFixed(d);
}

export function couplingText(p: NmrPeak): string {
  return p.couplings.length ? `J = ${p.couplings.map((c) => c.j.toFixed(1)).join(', ')} Hz` : '';
}
