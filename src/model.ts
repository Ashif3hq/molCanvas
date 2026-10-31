// Core document model for MolCanvas.
// Coordinates are in screen pixels (y grows downwards); one standard bond is BOND px long.

export const BOND = 40;

export type BondStereo = 'none' | 'wedge' | 'hash' | 'wavy';
export type ArrowKind = 'forward' | 'equilibrium' | 'retro' | 'resonance';
/** Non-bonding electrons drawn as dots: one unpaired (•), two unpaired (triplet ••) or two paired (singlet ••). */
export type Radical = 'doublet' | 'triplet' | 'singlet';

export interface Atom {
  id: number;
  el: string; // element symbol, or '*' for a polymer end / attachment point
  x: number;
  y: number;
  charge: number;
  iso?: number;
  radical?: Radical;
  sub?: string; // subscript shown next to a polymer bracket (default "n")
}

export interface Bond {
  id: number;
  a: number;
  b: number;
  order: 1 | 2 | 3;
  stereo: BondStereo;
}

export interface Arrow {
  id: number;
  x1: number;
  y1: number;
  x2: number;
  y2: number;
  kind: ArrowKind;
}

export interface TextItem {
  id: number;
  x: number;
  y: number;
  text: string;
}

/**
 * A repeat-unit bracket pair drawn around a block of atoms inside a chain, e.g. each block of a
 * block copolymer –[A]ₙ–[B]ₘ–. The two brackets cross the two bonds that leave the block.
 */
export interface BracketGroup {
  id: number;
  atoms: number[];
  sub: string; // count shown on the right-hand bracket ("n", "m", "50"…)
}

export interface Doc {
  atoms: Atom[];
  bonds: Bond[];
  arrows: Arrow[];
  texts: TextItem[];
  brackets?: BracketGroup[];
  nextId: number;
}

export interface Point {
  x: number;
  y: number;
}

export function emptyDoc(): Doc {
  return { atoms: [], bonds: [], arrows: [], texts: [], brackets: [], nextId: 1 };
}

export function cloneDoc(d: Doc): Doc {
  return JSON.parse(JSON.stringify(d));
}

export function newId(d: Doc): number {
  return d.nextId++;
}

export function addAtom(d: Doc, el: string, x: number, y: number, charge = 0): Atom {
  const a: Atom = { id: newId(d), el, x, y, charge };
  d.atoms.push(a);
  return a;
}

export function addBond(d: Doc, a: number, b: number, order: 1 | 2 | 3 = 1, stereo: BondStereo = 'none'): Bond | null {
  if (a === b) return null;
  const existing = bondBetween(d, a, b);
  if (existing) return existing;
  const bond: Bond = { id: newId(d), a, b, order, stereo };
  d.bonds.push(bond);
  return bond;
}

export function atomById(d: Doc, id: number): Atom | undefined {
  return d.atoms.find((a) => a.id === id);
}

export function bondBetween(d: Doc, a: number, b: number): Bond | undefined {
  return d.bonds.find((bd) => (bd.a === a && bd.b === b) || (bd.a === b && bd.b === a));
}

export function bondsOf(d: Doc, atomId: number): Bond[] {
  return d.bonds.filter((b) => b.a === atomId || b.b === atomId);
}

export function neighbors(d: Doc, atomId: number): Atom[] {
  const out: Atom[] = [];
  for (const b of d.bonds) {
    const other = b.a === atomId ? b.b : b.b === atomId ? b.a : -1;
    if (other >= 0) {
      const at = atomById(d, other);
      if (at) out.push(at);
    }
  }
  return out;
}

export function removeAtoms(d: Doc, ids: Set<number>) {
  d.atoms = d.atoms.filter((a) => !ids.has(a.id));
  d.bonds = d.bonds.filter((b) => !ids.has(b.a) && !ids.has(b.b));
  if (d.brackets) {
    for (const g of d.brackets) g.atoms = g.atoms.filter((id) => !ids.has(id));
    d.brackets = d.brackets.filter((g) => g.atoms.length);
  }
}

/** Remove any objects (atoms, bonds, arrows, texts) whose id is in the set. */
export function removeObjects(d: Doc, ids: Set<number>) {
  removeAtoms(d, ids);
  d.bonds = d.bonds.filter((b) => !ids.has(b.id));
  d.arrows = d.arrows.filter((a) => !ids.has(a.id));
  d.texts = d.texts.filter((t) => !ids.has(t.id));
  if (d.brackets) d.brackets = d.brackets.filter((g) => !ids.has(g.id));
}

/** Carbon atoms left with no bonds after an erase are removed, like ChemDraw does. */
export function pruneLoneCarbons(d: Doc, candidates: Iterable<number>) {
  const lone = new Set<number>();
  for (const id of candidates) {
    const a = atomById(d, id);
    if (a && (a.el === 'C' || a.el === '*') && a.charge === 0 && !a.radical && bondsOf(d, id).length === 0) lone.add(id);
  }
  if (lone.size) removeAtoms(d, lone);
}

/** Connected components, as lists of atom ids. */
export function fragments(d: Doc): number[][] {
  const adj = new Map<number, number[]>();
  for (const a of d.atoms) adj.set(a.id, []);
  for (const b of d.bonds) {
    adj.get(b.a)?.push(b.b);
    adj.get(b.b)?.push(b.a);
  }
  const seen = new Set<number>();
  const out: number[][] = [];
  for (const a of d.atoms) {
    if (seen.has(a.id)) continue;
    const comp: number[] = [];
    const stack = [a.id];
    seen.add(a.id);
    while (stack.length) {
      const cur = stack.pop()!;
      comp.push(cur);
      for (const n of adj.get(cur) ?? []) {
        if (!seen.has(n)) {
          seen.add(n);
          stack.push(n);
        }
      }
    }
    out.push(comp);
  }
  return out;
}

/** Copy the given atoms (and the bonds between them) into a new document. */
export function subDoc(d: Doc, atomIds: Set<number>, objIds?: Set<number>): Doc {
  const out = emptyDoc();
  out.atoms = d.atoms.filter((a) => atomIds.has(a.id)).map((a) => ({ ...a }));
  out.bonds = d.bonds.filter((b) => atomIds.has(b.a) && atomIds.has(b.b)).map((b) => ({ ...b }));
  if (objIds) {
    out.arrows = d.arrows.filter((a) => objIds.has(a.id)).map((a) => ({ ...a }));
    out.texts = d.texts.filter((t) => objIds.has(t.id)).map((t) => ({ ...t }));
  }
  out.brackets = (d.brackets ?? []).filter((g) => g.atoms.every((id) => atomIds.has(id))).map((g) => ({ ...g, atoms: [...g.atoms] }));
  out.nextId = d.nextId;
  return out;
}

/**
 * Merge a fragment document into `d`, translating it by (dx, dy).
 * Returns a map from fragment ids to the new ids in `d`.
 */
export function mergeInto(d: Doc, frag: Doc, dx = 0, dy = 0, atomMap = new Map<number, number>()): Map<number, number> {
  const map = new Map<number, number>(atomMap);
  for (const a of frag.atoms) {
    if (map.has(a.id)) continue;
    const na = addAtom(d, a.el, a.x + dx, a.y + dy, a.charge);
    if (a.iso) na.iso = a.iso;
    if (a.radical) na.radical = a.radical;
    if (a.sub) na.sub = a.sub;
    map.set(a.id, na.id);
  }
  for (const b of frag.bonds) {
    const a = map.get(b.a)!;
    const c = map.get(b.b)!;
    if (a === c) continue;
    const ex = bondBetween(d, a, c);
    if (ex) continue;
    const nb = addBond(d, a, c, b.order, b.stereo);
    if (nb) map.set(b.id, nb.id);
  }
  for (const ar of frag.arrows) {
    const id = newId(d);
    d.arrows.push({ ...ar, id, x1: ar.x1 + dx, y1: ar.y1 + dy, x2: ar.x2 + dx, y2: ar.y2 + dy });
    map.set(ar.id, id);
  }
  for (const t of frag.texts) {
    const id = newId(d);
    d.texts.push({ ...t, id, x: t.x + dx, y: t.y + dy });
    map.set(t.id, id);
  }
  for (const g of frag.brackets ?? []) {
    const atoms = g.atoms.map((id) => map.get(id)).filter((id): id is number => id !== undefined);
    if (!atoms.length) continue;
    (d.brackets ??= []).push({ id: newId(d), atoms, sub: g.sub });
  }
  return map;
}

export function bounds(d: Doc, pad = 0): { x: number; y: number; w: number; h: number } | null {
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  const take = (x: number, y: number) => {
    x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y);
  };
  for (const a of d.atoms) take(a.x, a.y);
  for (const a of d.arrows) { take(a.x1, a.y1); take(a.x2, a.y2); }
  for (const t of d.texts) { take(t.x, t.y - 12); take(t.x + t.text.length * 8, t.y + 4); }
  if (x0 === Infinity) return null;
  return { x: x0 - pad, y: y0 - pad, w: x1 - x0 + 2 * pad, h: y1 - y0 + 2 * pad };
}

export function centroid(atoms: Point[]): Point {
  if (!atoms.length) return { x: 0, y: 0 };
  let x = 0, y = 0;
  for (const a of atoms) { x += a.x; y += a.y; }
  return { x: x / atoms.length, y: y / atoms.length };
}

export const dist = (a: Point, b: Point) => Math.hypot(a.x - b.x, a.y - b.y);

/** Distance from p to segment ab. */
export function segDist(p: Point, a: Point, b: Point): number {
  const dx = b.x - a.x, dy = b.y - a.y;
  const len2 = dx * dx + dy * dy || 1;
  let t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2;
  t = Math.max(0, Math.min(1, t));
  return Math.hypot(p.x - (a.x + t * dx), p.y - (a.y + t * dy));
}

/**
 * Polymer ends that form a repeat unit: '*' atoms with one bond, in molecules that have exactly
 * two of them. A lone '*' is just an attachment point and is drawn as a plain asterisk.
 */
export function repeatUnitEnds(d: Doc): Set<number> {
  const out = new Set<number>();
  // a molecule with in-chain brackets keeps its * ends as plain asterisks outside the brackets
  const bracketed = new Set(validBrackets(d).flatMap((c) => c.group.atoms));
  const degree = new Map<number, number>();
  for (const b of d.bonds) { degree.set(b.a, (degree.get(b.a) ?? 0) + 1); degree.set(b.b, (degree.get(b.b) ?? 0) + 1); }
  const byId = new Map(d.atoms.map((a) => [a.id, a]));
  for (const comp of fragments(d)) {
    if (comp.some((id) => bracketed.has(id))) continue;
    const ends = comp.filter((id) => byId.get(id)!.el === '*' && degree.get(id) === 1);
    if (ends.length === 2) ends.forEach((id) => out.add(id));
  }
  return out;
}

/** Electrons a radical takes out of bonding: 1 for a doublet, 2 for a carbene/nitrene (singlet or triplet). */
export function radicalElectrons(r?: Radical): number {
  return r === 'doublet' ? 1 : r ? 2 : 0;
}

/** Unpaired electrons: 1 for a doublet, 2 for a triplet, 0 for a singlet (its two electrons are paired). */
export function unpairedElectrons(r?: Radical): number {
  return r === 'doublet' ? 1 : r === 'triplet' ? 2 : 0;
}

export interface Crossing { bond: Bond; inner: Atom; outer: Atom }

/** The bonds a bracket group crosses (one end inside the group, the other outside). */
export function bracketCrossings(d: Doc, g: BracketGroup): Crossing[] {
  const inside = new Set(g.atoms);
  const byId = new Map(d.atoms.map((a) => [a.id, a]));
  const out: Crossing[] = [];
  for (const b of d.bonds) {
    const ia = inside.has(b.a), ib = inside.has(b.b);
    if (ia === ib) continue;
    const inner = byId.get(ia ? b.a : b.b), outer = byId.get(ia ? b.b : b.a);
    if (inner && outer) out.push({ bond: b, inner, outer });
  }
  return out;
}

/** Bracket groups that can be drawn: a connected block of atoms left by exactly two bonds. */
export function validBrackets(d: Doc): { group: BracketGroup; crossings: Crossing[] }[] {
  return (d.brackets ?? []).map((group) => ({ group, crossings: bracketCrossings(d, group) }))
    .filter((c) => c.crossings.length === 2);
}

/**
 * Why a set of atoms cannot be a repeat unit, or null when it can: it must be one connected block
 * with exactly two bonds leading out of it (the bonds the brackets cross).
 */
export function bracketProblem(d: Doc, atomIds: Set<number>): string | null {
  if (!atomIds.size) return 'Drag a box around the atoms of one repeat unit.';
  const crossing = d.bonds.filter((b) => atomIds.has(b.a) !== atomIds.has(b.b)).length;
  const sub = subDoc(d, atomIds);
  if (fragments(sub).length > 1) return 'The box holds atoms that are not connected to each other. Box one repeat unit at a time.';
  if (crossing !== 2) {
    return crossing < 2
      ? 'A repeat unit needs two bonds leading out of the box, one on each side, for the brackets to cross. Draw the chain bonds (or * ends) on both sides first.'
      : `${crossing} bonds lead out of the box; a repeat unit needs exactly two. Make the box cover whole side groups.`;
  }
  return null;
}
