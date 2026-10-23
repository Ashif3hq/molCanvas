// Drawing geometry: where new bonds, rings, chains and templates go.
import { BondStereo, Atom, BOND, Doc, Point, atomById, bondBetween, bondsOf, emptyDoc, neighbors, dist } from './model';

const TAU = Math.PI * 2;
const MERGE_TOL = BOND * 0.25;

export const dirVec = (ang: number): Point => ({ x: Math.cos(ang), y: Math.sin(ang) });

function crowding(doc: Doc, p: Point, ignore: Set<number>): number {
  let s = 0;
  for (const a of doc.atoms) {
    if (ignore.has(a.id)) continue;
    const d = dist(a, p);
    if (d < BOND * 1.6) s += (BOND * 1.6 - d) / BOND;
  }
  return s;
}

/** Best angle (radians) for a new substituent on `atomId`. */
export function freeDirection(doc: Doc, atomId: number): number {
  const A = atomById(doc, atomId)!;
  const nb = neighbors(doc, atomId);
  if (nb.length === 0) return -Math.PI / 6;
  const angles = nb.map((n) => Math.atan2(n.y - A.y, n.x - A.x));
  if (nb.length === 1) {
    const th = angles[0];
    // sp centres (triple bond or cumulated double bonds) continue in a straight line
    const bs = bondsOf(doc, atomId);
    if (bs.some((b) => b.order === 3) || bs.filter((b) => b.order === 2).length >= 2) return th + Math.PI;
    const cands = [th + TAU / 3, th - TAU / 3];
    const N = nb[0];
    const grand = neighbors(doc, N.id).filter((g) => g.id !== A.id);
    const ignore = new Set([A.id]);
    let best = cands[0], bestScore = -Infinity;
    for (const c of cands) {
      const v = dirVec(c);
      const p = { x: A.x + v.x * BOND, y: A.y + v.y * BOND };
      let score = -crowding(doc, p, ignore) * 3;
      if (grand.length) {
        // trans zig-zag: the new bond should be parallel to grand -> N
        const g = grand[0];
        const gx = N.x - g.x, gy = N.y - g.y;
        const gl = Math.hypot(gx, gy) || 1;
        score += (v.x * gx + v.y * gy) / gl;
      } else {
        score += v.x * 0.5 + v.y * 0.01;
      }
      if (score > bestScore) { bestScore = score; best = c; }
    }
    return best;
  }
  const sorted = angles.map((a) => (a + TAU) % TAU).sort((a, b) => a - b);
  let bestGap = -1, bestAng = 0;
  for (let i = 0; i < sorted.length; i++) {
    const a = sorted[i];
    const b = i + 1 < sorted.length ? sorted[i + 1] : sorted[0] + TAU;
    if (b - a > bestGap) { bestGap = b - a; bestAng = (a + b) / 2; }
  }
  return bestAng;
}

export function snapAngle(ang: number, stepDeg = 30): number {
  const s = (stepDeg * Math.PI) / 180;
  return Math.round(ang / s) * s;
}

export function findAtomNear(doc: Doc, p: Point, tol = MERGE_TOL, exclude?: Set<number>): Atom | undefined {
  let best: Atom | undefined, bd = tol;
  for (const a of doc.atoms) {
    if (exclude?.has(a.id)) continue;
    const d = dist(a, p);
    if (d < bd) { bd = d; best = a; }
  }
  return best;
}

/** A fragment to merge: new atoms/bonds in `frag`, with some fragment atoms mapped onto existing atoms. */
export interface Placement {
  frag: Doc;
  map: Map<number, number>; // fragment atom id -> existing atom id
}

/** Map fragment atoms that land on top of existing atoms onto them. */
export function snapOntoExisting(doc: Doc, pl: Placement): Placement {
  const used = new Set(pl.map.values());
  for (const a of pl.frag.atoms) {
    if (pl.map.has(a.id) || a.el !== 'C') continue;
    const hit = findAtomNear(doc, a, MERGE_TOL, used);
    if (hit) { pl.map.set(a.id, hit.id); used.add(hit.id); }
  }
  return pl;
}

function ringVertices(center: Point, start: Point, n: number, dirSign: number): Point[] {
  const r = Math.hypot(start.x - center.x, start.y - center.y);
  const a0 = Math.atan2(start.y - center.y, start.x - center.x);
  const out: Point[] = [];
  for (let i = 0; i < n; i++) {
    const a = a0 + dirSign * (TAU * i) / n;
    out.push({ x: center.x + r * Math.cos(a), y: center.y + r * Math.sin(a) });
  }
  return out;
}

const ringRadius = (n: number) => BOND / (2 * Math.sin(Math.PI / n));
const ringApothem = (n: number) => BOND / (2 * Math.tan(Math.PI / n));

function buildRing(doc: Doc, pts: Point[], preset: Map<number, number>, aromatic: boolean): Placement {
  const frag = emptyDoc();
  frag.nextId = 1;
  const ids = pts.map((p) => {
    const id = frag.nextId++;
    frag.atoms.push({ id, el: 'C', x: p.x, y: p.y, charge: 0 });
    return id;
  });
  const map = new Map<number, number>();
  for (const [i, existing] of preset) map.set(ids[i], existing);
  const pl = snapOntoExisting(doc, { frag, map });
  const n = ids.length;
  const ringBonds: { a: number; b: number; existingOrder: number | null }[] = [];
  for (let i = 0; i < n; i++) {
    const a = ids[i], b = ids[(i + 1) % n];
    const ea = pl.map.get(a), eb = pl.map.get(b);
    const ex = ea !== undefined && eb !== undefined ? bondBetween(doc, ea, eb) : undefined;
    ringBonds.push({ a, b, existingOrder: ex ? ex.order : null });
  }
  const doubles = new Set<number>();
  if (aromatic) {
    // Atoms that already carry a double bond cannot take another one.
    const busy = new Set<number>();
    for (const id of ids) {
      const e = pl.map.get(id);
      if (e !== undefined && bondsOf(doc, e).some((b) => b.order >= 2)) busy.add(id);
    }
    const free = ringBonds.map((rb, i) => (rb.existingOrder === null ? i : -1)).filter((i) => i >= 0);
    let best: number[] = [];
    const total = 1 << free.length;
    for (let mask = 0; mask < total; mask++) {
      const chosen = free.filter((_, k) => mask & (1 << k));
      if (chosen.length <= best.length) continue;
      const usedAtoms = new Set<number>(busy);
      let ok = true;
      for (const i of chosen) {
        const { a, b } = ringBonds[i];
        if (usedAtoms.has(a) || usedAtoms.has(b)) { ok = false; break; }
        usedAtoms.add(a); usedAtoms.add(b);
      }
      if (ok) best = chosen;
    }
    best.forEach((i) => doubles.add(i));
  }
  ringBonds.forEach((rb, i) => {
    if (rb.existingOrder !== null) return;
    frag.bonds.push({ id: frag.nextId++, a: rb.a, b: rb.b, order: doubles.has(i) ? 2 : 1, stereo: 'none' });
  });
  return pl;
}

export function ringAtPoint(doc: Doc, p: Point, n: number, aromatic: boolean): Placement {
  const top = { x: p.x, y: p.y - ringRadius(n) };
  return buildRing(doc, ringVertices(p, top, n, 1), new Map(), aromatic);
}

export function ringOnAtom(doc: Doc, atomId: number, n: number, aromatic: boolean): Placement {
  const A = atomById(doc, atomId)!;
  const nb = neighbors(doc, atomId);
  let ang: number;
  if (nb.length === 0) ang = Math.PI / 2;
  else if (nb.length === 1) ang = Math.atan2(A.y - nb[0].y, A.x - nb[0].x);
  else ang = freeDirection(doc, atomId);
  const R = ringRadius(n);
  const c = { x: A.x + Math.cos(ang) * R, y: A.y + Math.sin(ang) * R };
  return buildRing(doc, ringVertices(c, A, n, 1), new Map([[0, atomId]]), aromatic);
}

export function ringOnBond(doc: Doc, bondId: number, n: number, aromatic: boolean, ringCenter?: Point, cursor?: Point): Placement {
  const b = doc.bonds.find((x) => x.id === bondId)!;
  const A = atomById(doc, b.a)!, B = atomById(doc, b.b)!;
  const mx = (A.x + B.x) / 2, my = (A.y + B.y) / 2;
  const len = dist(A, B) || 1;
  const px = -(B.y - A.y) / len, py = (B.x - A.x) / len;
  let side: number;
  if (ringCenter) {
    side = (ringCenter.x - mx) * px + (ringCenter.y - my) * py > 0 ? -1 : 1;
  } else {
    let score = 0;
    for (const nbr of [...neighbors(doc, A.id), ...neighbors(doc, B.id)]) {
      if (nbr.id === A.id || nbr.id === B.id) continue;
      score += Math.sign((nbr.x - mx) * px + (nbr.y - my) * py);
    }
    if (score === 0 && cursor) score = -Math.sign((cursor.x - mx) * px + (cursor.y - my) * py) || -1;
    side = score > 0 ? -1 : 1;
  }
  const ap = ringApothem(n) * (len / BOND);
  const c = { x: mx + px * ap * side, y: my + py * ap * side };
  // Walk from A to B around the centre.
  const cross = (A.x - c.x) * (B.y - c.y) - (A.y - c.y) * (B.x - c.x);
  const pts = ringVertices(c, A, n, cross > 0 ? 1 : -1);
  return buildRing(doc, pts, new Map([[0, A.id], [1, B.id]]), aromatic);
}

/** Zig-zag chain from `start` towards `to`. `startAtom` is reused as the first atom if given. */
export function chain(doc: Doc, start: Point, to: Point, startAtom?: number): Placement & { count: number } {
  const frag = emptyDoc();
  frag.nextId = 1;
  const map = new Map<number, number>();
  const ang = Math.atan2(to.y - start.y, to.x - start.x);
  const step = BOND * Math.cos(Math.PI / 6);
  const count = Math.max(1, Math.round(dist(start, to) / step));
  let sign = 1;
  if (startAtom !== undefined) {
    const nb = neighbors(doc, startAtom);
    if (nb.length) {
      const S = atomById(doc, startAtom)!;
      const p1 = { x: S.x + Math.cos(ang + Math.PI / 6) * BOND, y: S.y + Math.sin(ang + Math.PI / 6) * BOND };
      const p2 = { x: S.x + Math.cos(ang - Math.PI / 6) * BOND, y: S.y + Math.sin(ang - Math.PI / 6) * BOND };
      const d1 = Math.min(...nb.map((n) => dist(n, p1))), d2 = Math.min(...nb.map((n) => dist(n, p2)));
      sign = d1 >= d2 ? 1 : -1;
    }
  }
  let cur = { ...start };
  let prev = frag.nextId++;
  frag.atoms.push({ id: prev, el: 'C', x: cur.x, y: cur.y, charge: 0 });
  if (startAtom !== undefined) map.set(prev, startAtom);
  for (let i = 0; i < count; i++) {
    const a = ang + sign * (Math.PI / 6) * (i % 2 === 0 ? 1 : -1);
    cur = { x: cur.x + Math.cos(a) * BOND, y: cur.y + Math.sin(a) * BOND };
    const id = frag.nextId++;
    frag.atoms.push({ id, el: 'C', x: cur.x, y: cur.y, charge: 0 });
    frag.bonds.push({ id: frag.nextId++, a: prev, b: id, order: 1, stereo: 'none' });
    prev = id;
  }
  return { ...snapOntoExisting(doc, { frag, map }), count };
}

/** One bond from an atom (or point) to a new atom at `end`, snapping onto existing atoms. */
export function singleBond(doc: Doc, start: Point, end: Point, startAtom: number | undefined, el: string, order: 1 | 2 | 3, stereo: BondStereo): Placement {
  const frag = emptyDoc();
  frag.nextId = 1;
  const a = frag.nextId++, b = frag.nextId++;
  frag.atoms.push({ id: a, el: startAtom === undefined ? 'C' : 'C', x: start.x, y: start.y, charge: 0 });
  frag.atoms.push({ id: b, el, x: end.x, y: end.y, charge: 0 });
  frag.bonds.push({ id: frag.nextId++, a, b, order, stereo });
  const map = new Map<number, number>();
  if (startAtom !== undefined) map.set(a, startAtom);
  else {
    const s = findAtomNear(doc, start);
    if (s) map.set(a, s.id);
  }
  const exclude = new Set(map.values());
  const e = findAtomNear(doc, end, MERGE_TOL, exclude);
  if (e) map.set(b, e.id);
  return { frag, map };
}

export function rotatePoint(p: Point, c: Point, ang: number): Point {
  const cos = Math.cos(ang), sin = Math.sin(ang);
  const dx = p.x - c.x, dy = p.y - c.y;
  return { x: c.x + dx * cos - dy * sin, y: c.y + dx * sin + dy * cos };
}

function bboxCenter(frag: Doc): Point {
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const a of frag.atoms) { x0 = Math.min(x0, a.x); y0 = Math.min(y0, a.y); x1 = Math.max(x1, a.x); y1 = Math.max(y1, a.y); }
  return { x: (x0 + x1) / 2, y: (y0 + y1) / 2 };
}

/**
 * Place a preset fragment. Groups ('*' attachment point) attach to `atomId` if given,
 * pointing along `direction` (defaults to the atom's free direction).
 */
export function placeTemplate(doc: Doc, src: Doc, kind: 'molecule' | 'group' | 'polymer', p: Point, atomId?: number, direction?: number): Placement {
  const frag: Doc = JSON.parse(JSON.stringify(src));
  const map = new Map<number, number>();
  if (kind === 'group') {
    const star = frag.atoms.find((a) => a.el === '*');
    const bond = star && frag.bonds.find((b) => b.a === star.id || b.b === star.id);
    if (star && bond) {
      const first = frag.atoms.find((a) => a.id === (bond.a === star.id ? bond.b : bond.a))!;
      if (atomId !== undefined) {
        const A = atomById(doc, atomId)!;
        const want = direction ?? freeDirection(doc, atomId);
        const have = Math.atan2(first.y - star.y, first.x - star.x);
        const rot = want - have;
        for (const a of frag.atoms) {
          const r = rotatePoint(a, star, rot);
          a.x = r.x - star.x + A.x;
          a.y = r.y - star.y + A.y;
        }
        // Mirror the group if that makes it less crowded against the existing structure.
        const mirrored = frag.atoms.map((a) => mirrorAcross(a, A, want));
        const ignore = new Set([atomId]);
        const crowd = (pts: Point[]) => pts.reduce((s, q) => s + crowding(doc, q, ignore), 0);
        if (crowd(mirrored) + 0.01 < crowd(frag.atoms)) frag.atoms.forEach((a, i) => { a.x = mirrored[i].x; a.y = mirrored[i].y; });
        map.set(star.id, atomId);
        return { frag, map };
      }
      // On empty canvas: drop the attachment point.
      frag.atoms = frag.atoms.filter((a) => a.id !== star.id);
      frag.bonds = frag.bonds.filter((b) => b.id !== bond.id);
    }
  }
  const c = bboxCenter(frag);
  for (const a of frag.atoms) { a.x += p.x - c.x; a.y += p.y - c.y; }
  return { frag, map };
}

function mirrorAcross(p: Point, origin: Point, ang: number): Point {
  const ux = Math.cos(ang), uy = Math.sin(ang);
  const dx = p.x - origin.x, dy = p.y - origin.y;
  const along = dx * ux + dy * uy;
  const perp = -dx * uy + dy * ux;
  return { x: origin.x + along * ux + perp * uy, y: origin.y + along * uy - perp * ux };
}

/**
 * Lay out a repeat unit the way chemists draw it: the two polymer ends ('*') left and right with
 * horizontal end bonds, so the brackets are vertical. Side groups on an end carbon are turned out
 * of the way (on PMMA, PIB or PTFE they point straight up and down).
 */
export function orientPolymer(doc: Doc): Doc {
  const stars = doc.atoms.filter((a) => a.el === '*' && neighbors(doc, a.id).length === 1);
  if (stars.length !== 2) return doc;
  levelStars(doc, stars[0], stars[1]);
  const path = shortestPath(doc, stars[0].id, stars[1].id);
  if (!path) return doc;
  const onPath = new Set(path);
  const right = stars[0].x > stars[1].x ? stars[0] : stars[1];
  const angDist = (x: number, y: number) => Math.abs(Math.atan2(Math.sin(x - y), Math.cos(x - y)));
  for (const s of stars) {
    const X = neighbors(doc, s.id)[0];
    const B = neighbors(doc, X.id).find((n) => onPath.has(n.id) && n.id !== s.id);
    if (!B) continue;
    const outward = s === right ? 0 : Math.PI;
    const bAng = Math.atan2(B.y - X.y, B.x - X.x);
    if (angDist(outward, bAng) < Math.PI / 3) continue; // no room for a horizontal end bond
    s.x = X.x + Math.cos(outward) * BOND;
    s.y = X.y + Math.sin(outward) * BOND;
    const sides = neighbors(doc, X.id).filter((n) => n.id !== s.id && n.id !== B.id);
    const subtrees = sides.map((n) => subtree(doc, n.id, X.id));
    if (!sides.length || sides.length > 2 || subtrees.some((t) => t === null)) continue;
    const angOf = (n: Atom) => Math.atan2(n.y - X.y, n.x - X.x);
    const turn = (k: number, target: number) => {
      const rot = target - angOf(sides[k]);
      for (const id of subtrees[k]!) {
        const at = atomById(doc, id)!;
        Object.assign(at, rotatePoint(at, X, rot));
      }
    };
    if (sides.length === 2) {
      // backbone through, side groups straight up and down (keep which one is on top)
      const [top, bottom] = sides[0].y <= sides[1].y ? [0, 1] : [1, 0];
      turn(top, -Math.PI / 2);
      turn(bottom, Math.PI / 2);
    } else {
      const cur = angOf(sides[0]);
      if (angDist(cur, outward) < Math.PI / 3 || angDist(cur, bAng) < Math.PI / 3) {
        const mid = Math.atan2(Math.sin(outward) + Math.sin(bAng), Math.cos(outward) + Math.cos(bAng));
        const opts = [mid + Math.PI, mid];
        turn(0, angDist(opts[0], cur) < angDist(opts[1], cur) ? opts[0] : opts[1]);
      }
    }
  }
  return doc;
}

function levelStars(doc: Doc, s1: Atom, s2: Atom) {
  let ang = Math.atan2(s2.y - s1.y, s2.x - s1.x);
  if (Math.abs(ang) > Math.PI / 2) ang += Math.PI; // keep the current left/right sense where possible
  const c = { x: (s1.x + s2.x) / 2, y: (s1.y + s2.y) / 2 };
  for (const a of doc.atoms) Object.assign(a, rotatePoint(a, c, -ang));
}

function shortestPath(doc: Doc, from: number, to: number): number[] | null {
  const prev = new Map<number, number>([[from, -1]]);
  const queue = [from];
  while (queue.length) {
    const cur = queue.shift()!;
    if (cur === to) break;
    for (const n of neighbors(doc, cur)) if (!prev.has(n.id)) { prev.set(n.id, cur); queue.push(n.id); }
  }
  if (!prev.has(to)) return null;
  const path: number[] = [];
  for (let c = to; c !== -1; c = prev.get(c)!) path.push(c);
  return path;
}

/** Atoms reachable from `start` without passing through `block`; null if the branch loops back to it. */
function subtree(doc: Doc, start: number, block: number): number[] | null {
  const seen = new Set([start]);
  const stack = [start];
  while (stack.length) {
    const cur = stack.pop()!;
    for (const n of neighbors(doc, cur)) {
      if (n.id === block) { if (cur !== start) return null; continue; }
      if (!seen.has(n.id)) { seen.add(n.id); stack.push(n.id); }
    }
  }
  return [...seen];
}
