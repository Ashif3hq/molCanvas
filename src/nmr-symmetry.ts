// Stereo-aware NMR equivalence: which atoms (and which hydrogens) give the same signal.
//
// Two nuclei are isochronous when some symmetry of the time-averaged structure swaps them. Free rotation
// about single bonds and ring flips are assumed fast, amide C–N rotation slow. Every hydrogen is turned
// into a labelled atom (tritium stand-in) so the canonical ranking sees it, then:
//
//   • Canonical symmetry ranks that distinguish stereoheterotopic atoms (OpenChemLib's CIP-style ranking
//     with stereoheterotopicity, the same machinery RDKit's CanonicalRankAtoms(breakTies=False) uses)
//     give the orbits under proper symmetry: atoms with equal rank are homotopic.
//   • If the molecule is identical to its mirror image, canonical numbering of the molecule and of its
//     mirror image gives an improper symmetry; atoms it swaps are enantiotopic and stay equivalent.
//   • Atoms with the same constitution (stereo ignored) but different classes are diastereotopic.
//
// Stereocentres drawn without wedges are tried in every relative configuration (up to four of them), and two atoms are
// equivalent only when they are equivalent in all of them: a CH₂ between two undrawn stereocentres
// (2,4-pentanediol) is diastereotopic in the meso form, so it is reported as two signals.
//
// References: K. Mislow & M. Raban, Top. Stereochem. 1 (1967) 1 (topicity); W. T. Wipke & T. M. Dyott,
// J. Am. Chem. Soc. 96 (1974) 4834 (stereochemically unique naming); T. Sander, OpenChemLib Canonizer.
import { Canonizer, Molecule } from 'openchemlib';

type Mol = InstanceType<typeof Molecule>;

export type Topicity = 'homotopic' | 'enantiotopic' | 'diastereotopic' | 'constitutional';

export interface Equivalence {
  /** NMR class of each atom of the input molecule (by input index). Same class = same signal. */
  cls: number[];
  /** Constitution-only class of each input atom (stereo ignored). */
  constitution: number[];
  /** Per input heavy atom: NMR class of each of its ¹H (implicit or drawn), sorted. */
  hCls: number[][];
  /** Per input heavy atom: constitution class of its ¹H (all equal for H on one atom). */
  hConst: number[][];
  /** Per input heavy atom: proper-symmetry rank of each of its ¹H, in the same order as hCls. */
  hProper: number[][];
  /** Proper-symmetry rank of each input atom (homotopic atoms share it). */
  proper: number[];
  /** Input atoms that are stereocentres (drawn or not). */
  stereo: Set<number>;
  /** Stereocentres whose configuration was not drawn and was enumerated. */
  undrawn: number;
}

const TRITIUM = 3;
/** Undrawn stereocentres beyond this many keep one arbitrary configuration (2^(n-1) isomers are tried). */
const MAX_ENUMERATED = 4;

/** Dense 0-based class numbers from arbitrary keys, in order of first appearance. */
function dense(keys: (string | number)[]): number[] {
  const m = new Map<string | number, number>();
  return keys.map((k) => { if (!m.has(k)) m.set(k, m.size); return m.get(k)!; });
}

function symmetryRanks(m: Mol): number[] {
  return new Canonizer(m, { createSymmetryRank: true, neglectAnyStereoInformation: true }).getSymmetryRanks();
}

/**
 * Substitution test for one atom: mark it, then encode the marked molecule twice, with absolute stereo and
 * with every stereocentre made relative (racemic). Two atoms whose marked molecules are identical are
 * homotopic; mirror images (equal only when racemic) are enantiotopic; anything else is heterotopic.
 */
function substitutionIds(m: Mol, i: number): [string, string] {
  const t = m.getCompactCopy();
  t.setAtomCustomLabel(i, `${t.getAtomLabel(i)}*`);
  t.setAtomMass(i, t.getAtomicNo(i) === 1 ? 4 : t.getAtomMass(i) + 5);
  t.invalidateHelperArrays(Molecule.cHelperAll);
  const opts = { encodeAtomCustomLabels: true };
  const abs = new Canonizer(t, opts).getIDCode();
  t.ensureHelperArrays(Molecule.cHelperParities);
  for (let a = 0; a < t.getAllAtoms(); a++) if (t.getAtomParity(a) !== Molecule.cAtomParityNone) t.setAtomESR(a, Molecule.cESRTypeAnd, 0);
  const rac = new Canonizer(t, opts).getIDCode();
  return [abs, rac];
}

/**
 * Equivalence classes for a molecule as built for the predictor (2D coordinates, wedge bonds, hydrogens
 * implicit or drawn). Deuterium is kept as a different isotope, so CHD groups are handled too.
 */
export function equivalence(src: Mol): Equivalence {
  const m = src.getCompactCopy();
  const n0 = m.getAllAtoms();
  for (let i = 0; i < n0; i++) m.setAtomMapNo(i, i + 1, false);
  m.addImplicitHydrogens();
  const n = m.getAllAtoms();
  for (let i = 0; i < n; i++) if (m.getAtomicNo(i) === 1 && !m.getAtomMass(i) && !m.getAtomCharge(i)) m.setAtomMass(i, TRITIUM);
  m.ensureHelperArrays(Molecule.cHelperRings);
  const orig = (i: number) => m.getAtomMapNo(i) - 1; // -1 for hydrogens added here
  const isProton = (i: number) => m.getAtomicNo(i) === 1 && m.getAtomMass(i) === TRITIUM;

  // Slow amide rotation: the two N-substituents of a tertiary amide (DMF, N,N-diethylamides) sit cis and
  // trans to the carbonyl O. Labelling one of them breaks the symmetry the way the partial C=N bond does.
  const constitution = symmetryRanks(m);
  for (let i = 0; i < n; i++) {
    if (m.getAtomicNo(i) !== 7 || m.getAtomCharge(i) || m.isRingAtom(i) || m.getConnAtoms(i) !== 3) continue;
    const nb = Array.from({ length: 3 }, (_, k) => m.getConnAtom(i, k));
    const acyl = nb.filter((a) => isCarbonyl(m, a));
    if (acyl.length !== 1) continue;
    const rest = nb.filter((a) => a !== acyl[0]);
    if (rest.length === 2 && constitution[rest[0]] === constitution[rest[1]] && m.getAtomicNo(rest[0]) === 6) m.setAtomMass(rest[0], 13);
  }

  // Every sp³ atom that is not a drawn stereocentre gets a wedge to one of its ligands. Wedges on atoms that
  // are not stereocentres mean nothing, but once a substitution makes the atom stereogenic they fix its
  // configuration, so "replace H_a" and "replace H_b" give comparable molecules (a 2D stand-in for a 3D model).
  m.ensureHelperArrays(Molecule.cHelperParities);
  const stereoAtoms: number[] = [];
  const unknownWedge: number[] = [];
  const bridgeWedge: number[] = [];
  const ringBonds = (a: number) => { let r = 0; for (let k = 0; k < m.getConnAtoms(a); k++) if (m.isRingBond(m.getConnBond(a, k))) r++; return r; };
  const hasStereoBond = (a: number) => {
    for (let k = 0; k < m.getAllConnAtoms(a); k++) {
      const b = m.getConnBond(a, k);
      if (m.isStereoBond(b) && m.getBondAtom(0, b) === a) return true;
    }
    return false;
  };
  const addWedge = (a: number): number => {
    // Prefer a hydrogen, then a non-ring single bond, then any plain single bond.
    const cands: number[] = [];
    for (let k = 0; k < m.getAllConnAtoms(a); k++) {
      const b = m.getConnBond(a, k), o = m.getConnAtom(a, k);
      if (m.getBondOrder(b) !== 1 || m.getBondType(b) !== Molecule.cBondTypeSingle) continue;
      const score = m.getAtomicNo(o) === 1 ? 0 : !m.isRingBond(b) ? 1 : 2;
      cands[score] ??= b;
    }
    const pick = cands.find((b) => b !== undefined) ?? -1;
    if (pick < 0) return -1;
    if (m.getBondAtom(0, pick) !== a) { const other = m.getBondAtom(0, pick); m.setBondAtom(1, pick, other); m.setBondAtom(0, pick, a); }
    m.setBondType(pick, Molecule.cBondTypeUp);
    return pick;
  };
  for (let i = 0; i < n; i++) {
    if (m.getAtomicNo(i) === 1) continue;
    if (m.isAtomStereoCenter(i)) {
      stereoAtoms.push(i);
      if (m.getAtomParity(i) === Molecule.cAtomParityUnknown) { const b = addWedge(i); if (b >= 0) unknownWedge.push(b); }
      continue;
    }
    if (m.getAllConnAtoms(i) !== 4 || hasStereoBond(i) || m.getAtomPi(i) !== 0) continue;
    // Bridgeheads (three or more ring bonds) that are not stereocentres (norbornane, adamantane) have their
    // configuration fixed by the cage; their helper wedges are tried both ways below, and the most symmetric
    // result is kept, since an arbitrary choice can describe a geometry that cannot exist.
    const b = addWedge(i);
    if (b >= 0 && ringBonds(i) >= 3) bridgeWedge.push(b);
  }

  // Only atoms that share their constitution with another atom need the (slow) substitution test.
  const size = new Map<number, number>();
  for (let i = 0; i < n; i++) size.set(constitution[i], (size.get(constitution[i]) ?? 0) + 1);
  const unique = (i: number) => size.get(constitution[i]) === 1;
  // Hydrogens of CH, CH₃ (fast rotation), OH and NH share the class of their atom; only the two H of an
  // XH₂ group can differ from each other, so only those are tested one by one.
  const protons = new Array<number>(n).fill(0);
  for (let i = 0; i < n; i++) if (isProton(i)) protons[m.getConnAtom(i, 0)]++;
  const follows = (i: number) => isProton(i) && protons[m.getConnAtom(i, 0)] !== 2;
  const setWedges = (bonds: number[], bits: number, count: number) => bonds.forEach((b, j) => {
    const up = j === 0 || j >= count ? true : ((bits >> (j - 1)) & 1) === 1;
    m.setBondType(b, up ? Molecule.cBondTypeUp : Molecule.cBondTypeDown);
  });

  // Undrawn stereocentres: every relative configuration (mirror images are covered by the racemic test);
  // atoms are equivalent only if they are equivalent in all of them.
  const k = Math.min(unknownWedge.length, MAX_ENUMERATED);
  const classes = (): { cls: number[]; proper: number[] } => {
    const keys: string[][] = Array.from({ length: n }, () => []);
    const properKeys: string[] = new Array(n).fill('');
    for (let c = 0; c < (k ? 1 << (k - 1) : 1); c++) {
      setWedges(unknownWedge, c, k);
      m.invalidateHelperArrays(Molecule.cHelperAll);
      for (let i = 0; i < n; i++) {
        if (follows(i)) continue;
        if (unique(i)) { keys[i].push(`u${constitution[i]}`); properKeys[i] = `u${constitution[i]}`; continue; }
        const [abs, rac] = substitutionIds(m, i);
        keys[i].push(rac);
        if (c === 0) properKeys[i] = abs;
      }
    }
    for (let i = 0; i < n; i++) if (follows(i)) { const a = m.getConnAtom(i, 0); keys[i] = keys[a].map((x) => `H on ${x}`); properKeys[i] = `H on ${properKeys[a]}`; }
    return { cls: dense(keys.map((kk) => kk.join(' '))), proper: dense(properKeys) };
  };
  let best = { cls: [] as number[], proper: [] as number[] };
  const kb = Math.min(bridgeWedge.length, MAX_ENUMERATED);
  for (let c = 0; c < (kb ? 1 << (kb - 1) : 1); c++) {
    setWedges(bridgeWedge, c, kb);
    const r = classes();
    if (!best.cls.length || new Set(r.cls).size < new Set(best.cls).size) best = r;
  }
  const { cls, proper } = best;

  const out: Equivalence = {
    cls: new Array(n0).fill(-1), constitution: new Array(n0).fill(-1), proper: new Array(n0).fill(-1),
    hCls: Array.from({ length: n0 }, () => []), hConst: Array.from({ length: n0 }, () => []), hProper: Array.from({ length: n0 }, () => []),
    stereo: new Set(), undrawn: unknownWedge.length,
  };
  for (let i = 0; i < n; i++) {
    const o = orig(i);
    if (o >= 0) { out.cls[o] = cls[i]; out.constitution[o] = constitution[i]; out.proper[o] = proper[i]; }
    if (isProton(i)) {
      const a = orig(m.getConnAtom(i, 0));
      if (a < 0) continue;
      out.hCls[a].push(cls[i]); out.hConst[a].push(constitution[i]); out.hProper[a].push(proper[i]);
    }
  }
  for (const a of stereoAtoms) { const o = orig(a); if (o >= 0) out.stereo.add(o); }
  return out;
}

function isCarbonyl(m: Mol, a: number): boolean {
  if (m.getAtomicNo(a) !== 6) return false;
  for (let k = 0; k < m.getConnAtoms(a); k++) {
    const b = m.getConnBond(a, k);
    if (m.getBondOrder(b) === 2 && m.getAtomicNo(m.getConnAtom(a, k)) === 8) return true;
  }
  return false;
}

/** Topicity of two atoms (or two hydrogens) from their classes. */
export function topicity(sameConstitution: boolean, sameClass: boolean, sameProper: boolean): Topicity {
  if (!sameConstitution) return 'constitutional';
  if (!sameClass) return 'diastereotopic';
  return sameProper ? 'homotopic' : 'enantiotopic';
}
