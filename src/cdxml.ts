// ChemDraw CDXML import and export.
// CDXML is ChemDraw's XML format: coordinates in points (1/72 inch, y down), atoms are <n> nodes and
// bonds are <b> elements inside <fragment>s. Abbreviations (OMe, Boc…) are nodes that carry their
// expanded structure as a nested fragment. Polymer repeat units are <bracketedgroup>s with SRU usage.
import { Molecule } from 'openchemlib';
import { analyze, toOCL } from './chem';
import { ArrowKind, Atom, BOND, Bond, BondStereo, Doc, emptyDoc, fragments, Point, Radical, repeatUnitEnds, TextItem } from './model';
import { bracketMarks, FONT, labelInfo, textBox, textWidth } from './render';

// ---------- a small XML reader (CDXML only needs elements, attributes and text) ----------

export interface XNode {
  tag: string;
  attrs: Record<string, string>;
  children: XNode[];
  text: string;
}

const ENTITIES: Record<string, string> = { amp: '&', lt: '<', gt: '>', quot: '"', apos: "'" };
const decode = (s: string) => s.replace(/&(#x[0-9a-f]+|#\d+|\w+);/gi, (m, e: string) => {
  if (e[0] === '#') return String.fromCodePoint(e[1] === 'x' || e[1] === 'X' ? parseInt(e.slice(2), 16) : parseInt(e.slice(1), 10));
  return ENTITIES[e] ?? m;
});

export function parseXml(src: string): XNode {
  const root: XNode = { tag: '#root', attrs: {}, children: [], text: '' };
  const stack = [root];
  const token = /<!--[\s\S]*?-->|<\?[\s\S]*?\?>|<!\[CDATA\[([\s\S]*?)\]\]>|<![^>]*>|<\/\s*([\w:.-]+)\s*>|<([\w:.-]+)((?:\s+[\w:.-]+\s*=\s*(?:"[^"]*"|'[^']*'))*)\s*(\/?)>|([^<]+)|(<)/g;
  const attr = /([\w:.-]+)\s*=\s*(?:"([^"]*)"|'([^']*)')/g;
  let m: RegExpExecArray | null;
  while ((m = token.exec(src))) {
    const top = stack[stack.length - 1];
    if (m[1] !== undefined) top.text += m[1];
    else if (m[2]) {
      if (top.tag !== m[2]) throw new Error(`Malformed XML: </${m[2]}> does not close <${top.tag}>`);
      stack.pop();
    } else if (m[3]) {
      const node: XNode = { tag: m[3], attrs: {}, children: [], text: '' };
      let a: RegExpExecArray | null;
      attr.lastIndex = 0;
      while ((a = attr.exec(m[4]))) node.attrs[a[1]] = decode(a[2] ?? a[3]);
      top.children.push(node);
      if (!m[5]) stack.push(node);
    } else if (m[6] !== undefined) top.text += decode(m[6]);
    else if (m[7]) throw new Error('Malformed XML: stray “<”');
  }
  if (stack.length > 1) throw new Error(`Malformed XML: <${stack[stack.length - 1].tag}> is never closed`);
  return root;
}

const num = (s?: string) => (s ?? '').split(/\s+/).filter(Boolean).map(Number);
const textOf = (t: XNode) => t.children.filter((c) => c.tag === 's').map((s) => s.text).join('') || t.text;

// ---------- import ----------

export interface CdxmlImport {
  doc: Doc;
  /** Plain-language notes about anything that could not be carried over exactly. */
  notes: string[];
}

interface RawNode { id: string; type: string; el: string; p?: Point; charge: number; iso?: number; radical?: Radical; numH?: number; label?: string; labelP?: Point; extNum?: number; inner?: string[] }
interface RawBond { id: string; b: string; e: string; order: string; display: string; beginExt?: number; endExt?: number }


export function fromCdxml(text: string): CdxmlImport {
  const root = parseXml(text).children.find((c) => c.tag === 'CDXML');
  if (!root) throw new Error('not a ChemDraw CDXML file');
  const nodes = new Map<string, RawNode>();
  const bonds: RawBond[] = [];
  const texts: { p: Point; text: string; just: string }[] = [];
  const arrows: { tail: Point; head: Point; kind: ArrowKind }[] = [];
  const pluses: Point[] = [];
  const groups: XNode[] = [];
  const skipped = new Map<string, number>();
  const skip = (what: string) => skipped.set(what, (skipped.get(what) ?? 0) + 1);
  const notes: string[] = [];
  const dropped = new Set<string>();

  const point = (s: string | undefined): Point | undefined => {
    const v = num(s);
    return v.length >= 2 && v.every(Number.isFinite) ? { x: v[0], y: v[1] } : undefined;
  };

  const readText = (t: XNode) => {
    const p = point(t.attrs.p);
    const s = textOf(t).replace(/\r/g, '\n');
    if (p && s.trim()) texts.push({ p, text: s, just: t.attrs.Justification ?? t.attrs.CaptionJustification ?? 'Left' });
  };

  const readGraphic = (g: XNode) => {
    if (g.attrs.SupersededBy) return; // an <arrow> with the same meaning follows
    const type = g.attrs.GraphicType;
    const bb = num(g.attrs.BoundingBox);
    if (type === 'Line') {
      const at = g.attrs.ArrowType ?? 'NoHead';
      const kind: ArrowKind | null = at === 'Equilibrium' ? 'equilibrium' : at === 'RetroSynthetic' ? 'retro' : at === 'Resonance' ? 'resonance'
        : /FullHead|HalfHead|Hollow/.test(at) ? 'forward' : null;
      if (kind && bb.length === 4) arrows.push({ head: { x: bb[0], y: bb[1] }, tail: { x: bb[2], y: bb[3] }, kind });
      else skip(kind ? 'arrows' : 'lines');
    } else if (type === 'Symbol' && (g.attrs.SymbolType === 'Plus') && bb.length >= 2) {
      pluses.push({ x: bb[0], y: bb[1] });
    } else if (type === 'Bracket') {
      // the bracket's meaning is read from <bracketedgroup>; MolCanvas draws its own brackets
    } else if (type === 'Arc' && g.attrs.ArrowType && g.attrs.ArrowType !== 'NoHead') {
      skip('curved arrows');
    } else if (type === 'Symbol' && g.children.some((c) => c.tag === 'represent' && c.attrs.attribute === 'Radical')) {
      // ChemDraw's radical dot, already carried by the node's Radical attribute
    } else if (type === 'Symbol' && g.attrs.SymbolType === 'LonePair') {
      // lone pairs follow from the structure; MolCanvas draws them with “Show lone pairs”
    } else if (type === 'Symbol' && /Electron|Radical/.test(g.attrs.SymbolType ?? '')) {
      skip('free-standing electron dots');
    } else if (type !== 'Symbol') {
      skip('shapes');
    }
  };

  const readArrow = (a: XNode) => {
    const head = point(a.attrs.Head3D), tail = point(a.attrs.Tail3D);
    if (!head || !tail) return;
    if (a.attrs.AngularSize && Math.abs(parseFloat(a.attrs.AngularSize)) > 0.5) { skip('curved arrows'); return; }
    const h = a.attrs.ArrowheadHead ?? 'Unspecified', t = a.attrs.ArrowheadTail ?? 'Unspecified';
    const type = a.attrs.ArrowheadType ?? 'Solid';
    let kind: ArrowKind;
    if (/Half/.test(h) && /Half/.test(t)) kind = 'equilibrium';
    else if (h === 'Full' && t === 'Full') kind = 'resonance';
    else if (type === 'Angle' && a.attrs.ArrowShaftSpacing) kind = 'retro';
    else if (/Full|Half/.test(h)) kind = 'forward';
    else if (/Full|Half/.test(t)) { arrows.push({ head: tail, tail: head, kind: 'forward' }); return; }
    else { skip('lines'); return; }
    if (a.attrs.NoGo) skip('“no reaction” crosses on arrows');
    arrows.push({ head, tail, kind });
  };

  const readNode = (n: XNode, inheritP?: Point): RawNode => {
    const type = n.attrs.NodeType ?? 'Element';
    const p = point(n.attrs.p) ?? inheritP;
    const t = n.children.find((c) => c.tag === 't');
    const node: RawNode = {
      id: n.attrs.id, type, el: 'C', p, charge: parseInt(n.attrs.Charge ?? '0', 10) || 0,
      label: t ? textOf(t) : n.attrs.GenericNickname, labelP: t ? point(t.attrs.p) : undefined,
    };
    if (type === 'Element') {
      const z = parseInt(n.attrs.Element ?? '6', 10);
      const sym = Molecule.cAtomLabel[z];
      node.el = z > 0 && sym ? sym : '*';
      if (n.attrs.Isotope) node.iso = parseInt(n.attrs.Isotope, 10);
      if (n.attrs.NumHydrogens !== undefined) node.numH = parseInt(n.attrs.NumHydrogens, 10);
      const rad = n.attrs.Radical;
      if (rad === 'Doublet' || rad === 'Triplet' || rad === 'Singlet') node.radical = rad.toLowerCase() as Radical;
      else if (rad && rad !== 'None') skip('radicals');
    } else if (type === 'ExternalConnectionPoint') {
      node.extNum = parseInt(n.attrs.ExternalConnectionNum ?? '0', 10) || undefined;
    } else if (type === 'Nickname' || type === 'Fragment') {
      const inner = n.children.find((c) => c.tag === 'fragment');
      const ids = inner ? readFragment(inner, p) : [];
      // without an external connection point the label cannot be attached, so it stays a generic label
      if (ids.some((i) => nodes.get(i)!.type === 'ExternalConnectionPoint')) node.inner = ids;
      else { node.type = 'Unspecified'; ids.forEach((i) => dropped.add(i)); }
    }
    nodes.set(node.id, node);
    return node;
  };

  /** Reads a fragment's nodes and bonds; returns its node ids. */
  function readFragment(f: XNode, inheritP?: Point): string[] {
    const ids: string[] = [];
    for (const c of f.children) {
      if (c.tag === 'n') ids.push(readNode(c, inheritP).id);
      else if (c.tag === 'b') {
        bonds.push({
          id: c.attrs.id, b: c.attrs.B, e: c.attrs.E, order: (c.attrs.Order ?? '1').trim(), display: c.attrs.Display ?? 'Solid',
          beginExt: c.attrs.BeginExternalNum ? parseInt(c.attrs.BeginExternalNum, 10) : undefined,
          endExt: c.attrs.EndExternalNum ? parseInt(c.attrs.EndExternalNum, 10) : undefined,
        });
      } else if (c.tag === 'graphic') readGraphic(c);
      else if (c.tag === 't') readText(c);
      else if (c.tag === 'arrow') readArrow(c);
      else if (c.tag === 'bracketedgroup') groups.push(c);
      else if (c.tag === 'fragment') ids.push(...readFragment(c, inheritP));
    }
    return ids;
  }

  const walk = (x: XNode) => {
    for (const c of x.children) {
      if (c.tag === 'fragment') readFragment(c);
      else if (c.tag === 't') readText(c);
      else if (c.tag === 'graphic') readGraphic(c);
      else if (c.tag === 'arrow') readArrow(c);
      else if (c.tag === 'bracketedgroup') groups.push(c);
      else if (c.tag === 'page' || c.tag === 'group') walk(c);
      else if (c.tag === 'curve' || c.tag === 'table' || c.tag === 'embeddedobject' || c.tag === 'altgroup' || c.tag === 'spectrum' || c.tag === 'tlcplate') skip('shapes');
    }
  };
  walk(root);

  // Abbreviations: a bond to the label really goes to the atom next to its external connection point.
  const bondsOfNode = new Map<string, RawBond[]>();
  for (const b of bonds) for (const end of [b.b, b.e]) { if (!bondsOfNode.has(end)) bondsOfNode.set(end, []); bondsOfNode.get(end)!.push(b); }
  const resolve = (id: string, ext: number | undefined, used: Map<string, number>): string | null => {
    for (let guard = 0; guard < 50; guard++) {
      const n = nodes.get(id);
      if (!n || !n.inner) return n ? id : null;
      const ecps = n.inner.map((i) => nodes.get(i)!).filter((x) => x.type === 'ExternalConnectionPoint')
        .sort((a, b) => (a.extNum ?? 0) - (b.extNum ?? 0));
      let ecp = ext !== undefined ? ecps.find((e) => e.extNum === ext) : undefined;
      if (!ecp) { const k = used.get(id) ?? 0; ecp = ecps[k]; used.set(id, k + 1); }
      if (!ecp) return null;
      const link = bondsOfNode.get(ecp.id)?.[0];
      if (!link) return null;
      id = link.b === ecp.id ? link.e : link.b;
      ext = undefined;
    }
    return null;
  };

  const doc = emptyDoc();
  const atomOf = new Map<string, Atom>();
  const generic = new Set<string>();
  for (const n of nodes.values()) {
    if (n.type === 'ExternalConnectionPoint' || n.inner || dropped.has(n.id)) continue;
    if (!n.p) continue;
    let el = n.el;
    if (n.type !== 'Element' && !bondsOfNode.has(n.id)) {
      // a free-standing label such as “THF” or “−78 °C” is text, not an atom
      if (n.label?.trim()) texts.push({ p: n.labelP ?? n.p, text: n.label, just: 'Left' });
      continue;
    }
    if (n.type !== 'Element') {
      el = '*';
      generic.add(n.label?.trim() || '?');
    }
    const a: Atom = { id: doc.nextId++, el, x: n.p.x, y: n.p.y, charge: el === '*' ? 0 : n.charge };
    if (n.iso) a.iso = n.iso;
    if (n.radical && el !== '*') a.radical = n.radical;
    doc.atoms.push(a);
    atomOf.set(n.id, a);
  }

  const bondOf = new Map<string, Bond>();
  const aromatic = new Set<number>();
  const used = new Map<string, number>();
  let datives = 0;
  for (const rb of bonds) {
    const nb = nodes.get(rb.b), ne = nodes.get(rb.e);
    if (!nb || !ne || nb.type === 'ExternalConnectionPoint' || ne.type === 'ExternalConnectionPoint') continue;
    const ib = resolve(rb.b, rb.beginExt, used), ie = resolve(rb.e, rb.endExt, used);
    const A = ib ? atomOf.get(ib) : undefined, B = ie ? atomOf.get(ie) : undefined;
    if (!A || !B || A === B) continue;
    if (doc.bonds.some((x) => (x.a === A.id && x.b === B.id) || (x.a === B.id && x.b === A.id))) continue;
    const ord = rb.order.split(/\s+/)[0];
    let order: 1 | 2 | 3;
    if (ord === '1' || ord === '2' || ord === '3') order = parseInt(ord, 10) as 1 | 2 | 3;
    else if (ord === '1.5') order = 1;
    else if (ord === 'dative') { order = 1; A.charge += 1; B.charge -= 1; datives++; }
    else { skip(ord === 'hydrogen' ? 'hydrogen bonds' : ord === 'ionic' ? 'ionic bonds' : 'bonds of unsupported type'); continue; }
    let a = A.id, b = B.id, stereo: BondStereo = 'none';
    if (order === 1) {
      const d = rb.display;
      if (d === 'WedgeBegin' || d === 'WedgeEnd') stereo = 'wedge';
      else if (d === 'WedgedHashBegin' || d === 'WedgedHashEnd') stereo = 'hash';
      else if (d === 'Wavy') stereo = 'wavy';
      if (d.endsWith('End')) { a = B.id; b = A.id; } // the narrow end is at E
    } else if (order === 2 && rb.display === 'Wavy') stereo = 'wavy'; // E/Z unknown
    const bond: Bond = { id: doc.nextId++, a, b, order, stereo };
    if (ord === '1.5') aromatic.add(bond.id);
    doc.bonds.push(bond);
    bondOf.set(rb.id, bond);
  }
  if (aromatic.size && !kekulize(doc, aromatic)) notes.push('Some aromatic bonds could not be assigned alternating single and double bonds and were left single.');

  // Polymer repeat units: the atoms outside the brackets become polymer ends (*).
  for (const g of groups) {
    const usage = g.attrs.BracketUsage ?? '';
    if (usage !== 'SRU') { skip(usage === 'MultipleGroup' ? 'multiple-group brackets' : 'brackets'); continue; }
    const crossing = g.children.flatMap((c) => c.children).filter((c) => c.tag === 'crossingbond');
    const ends: Atom[] = [];
    for (const cb of crossing) {
      const bond = bondOf.get(cb.attrs.BondID);
      const inner = cb.attrs.InnerAtomID ? atomOf.get(resolve(cb.attrs.InnerAtomID, undefined, new Map()) ?? '') : undefined;
      if (!bond || !inner) continue;
      const outerId = bond.a === inner.id ? bond.b : bond.a;
      const outer = doc.atoms.find((x) => x.id === outerId)!;
      const degree = doc.bonds.filter((x) => x.a === outer.id || x.b === outer.id).length;
      if (degree === 1) ends.push(outer);
    }
    if (ends.length === 2 && crossing.length === 2) {
      for (const e of ends) { e.el = '*'; e.charge = 0; delete e.iso; e.sub = g.attrs.SRULabel || 'n'; }
      for (const b of doc.bonds) if (ends.some((e) => e.id === b.a || e.id === b.b)) { b.order = 1; b.stereo = 'none'; }
    } else if (crossing.length === 2) {
      // a block inside a longer chain (e.g. one block of a copolymer) keeps its brackets and count
      const atoms = (g.attrs.BracketedObjectIDs ?? '').split(/\s+/).map((id) => atomOf.get(id)?.id).filter((id): id is number => id !== undefined);
      if (atoms.length) (doc.brackets ??= []).push({ id: doc.nextId++, atoms, sub: g.attrs.SRULabel || 'n' });
    } else {
      notes.push('A repeat-unit bracket that crosses more or fewer than two bonds was not converted.');
    }
  }

  // Scale: ChemDraw's bond length (often 14.4 or 30 pt) becomes MolCanvas's standard bond.
  const lens = doc.bonds.map((b) => {
    const A = doc.atoms.find((x) => x.id === b.a)!, B = doc.atoms.find((x) => x.id === b.b)!;
    return Math.hypot(A.x - B.x, A.y - B.y);
  }).filter((l) => l > 0.01).sort((p, q) => p - q);
  const ref = lens.length ? lens[Math.floor(lens.length / 2)] : parseFloat(root.attrs.BondLength ?? '') || 14.4;
  const k = BOND / ref;
  for (const a of doc.atoms) { a.x *= k; a.y *= k; }
  for (const ar of arrows) {
    doc.arrows.push({ id: doc.nextId++, x1: ar.tail.x * k, y1: ar.tail.y * k, x2: ar.head.x * k, y2: ar.head.y * k, kind: ar.kind });
  }
  const lineH = FONT * 1.15;
  for (const t of texts) {
    t.text.split('\n').forEach((line, i) => {
      if (!line.trim()) return;
      const item: TextItem = { id: doc.nextId++, x: t.p.x * k, y: t.p.y * k + i * lineH, text: line.trim() };
      const w = textBox(item).w - 4;
      if (t.just === 'Center') item.x -= w / 2;
      else if (t.just === 'Right') item.x -= w;
      doc.texts.push(item);
    });
  }
  for (const p of pluses) {
    const q = { x: p.x * k, y: p.y * k };
    // a + next to an atom is ChemDraw's charge symbol, already carried by the atom's charge
    if (doc.atoms.some((a) => Math.hypot(a.x - q.x, a.y - q.y) < BOND * 0.8)) continue;
    doc.texts.push({ id: doc.nextId++, x: q.x - textWidth('+', FONT * 0.9) * 0.45, y: q.y + FONT * 0.32, text: '+' });
  }

  if (!doc.atoms.length && !doc.arrows.length && !doc.texts.length) throw new Error('no structures, arrows or text found');

  // Hydrogen counts ChemDraw stored that differ from the standard valence (radicals, unusual valences).
  const an = analyze(doc);
  let hMismatch = 0;
  for (const [id, a] of atomOf) {
    const n = nodes.get(id)!;
    if (n.numH !== undefined && a.el !== '*' && (an.hCount.get(a.id) ?? 0) !== n.numH) hMismatch++;
  }
  if (datives) notes.push(`${datives} dative bond${datives > 1 ? 's were' : ' was'} drawn as single bonds with formal charges (+ on the donor, − on the acceptor).`);
  if (generic.size) notes.push(`Generic labels (${[...generic].slice(0, 4).join(', ')}) became attachment points (*), since they do not stand for one structure.`);
  if (hMismatch) notes.push(`${hMismatch} atom${hMismatch > 1 ? 's have' : ' has'} a hydrogen count in the file that differs from its standard valence; MolCanvas uses the standard valence.`);
  if (skipped.size) notes.push(`Not supported and left out: ${[...skipped].map(([w, c]) => `${c} ${w}`).join(', ')}.`);
  return { doc, notes };
}

/** Turns aromatic (1.5) bonds into alternating single and double bonds. False if that is impossible. */
function kekulize(doc: Doc, aromatic: Set<number>): boolean {
  const { mol, bondIndex } = toOCL(doc);
  for (const id of aromatic) mol.setBondType(bondIndex.get(id)!, Molecule.cBondTypeDelocalized);
  for (let i = 0; i < mol.getAllAtoms(); i++) mol.setAtomMapNo(i, i + 1, false);
  let k: InstanceType<typeof Molecule>;
  try {
    k = Molecule.fromSmiles(mol.toIsomericSmiles({ kekulizedOutput: true, includeMapping: true }));
  } catch { return false; }
  const byMap = new Map<number, number>();
  for (let i = 0; i < k.getAllAtoms(); i++) byMap.set(k.getAtomMapNo(i), i);
  let ok = true;
  for (const b of doc.bonds) {
    if (!aromatic.has(b.id)) continue;
    const bi = bondIndex.get(b.id)!;
    const ka = byMap.get(mol.getBondAtom(0, bi) + 1), kb = byMap.get(mol.getBondAtom(1, bi) + 1);
    const kbond = ka !== undefined && kb !== undefined ? k.getBond(ka, kb) : -1;
    const order = kbond >= 0 ? k.getBondOrder(kbond) : 0;
    if (order === 1 || order === 2) b.order = order;
    else ok = false;
  }
  return ok;
}

// ---------- export ----------

const PT = 14.4 / BOND; // ACS Document 1996 style: 14.4 pt bonds, 10 pt Arial labels
const LABEL_SIZE = 10;
const xmlEsc = (s: string) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const r2 = (n: number) => (Math.round(n * 100) / 100).toFixed(2);

const ARROW_TYPE: Record<ArrowKind, string> = { forward: 'FullHead', equilibrium: 'Equilibrium', retro: 'RetroSynthetic', resonance: 'Resonance' };

/** A CDXML document ChemDraw opens with ACS 1996 settings. */
export function toCdxml(doc: Doc): string {
  let nextId = Math.max(doc.nextId, 1) + 1;
  const newId = () => nextId++;
  const an = analyze(doc);
  const byId = new Map(doc.atoms.map((a) => [a.id, a]));
  const nbrs = new Map<number, Atom[]>(doc.atoms.map((a) => [a.id, []]));
  for (const b of doc.bonds) { nbrs.get(b.a)?.push(byId.get(b.b)!); nbrs.get(b.b)?.push(byId.get(b.a)!); }
  const ends = repeatUnitEnds(doc);

  // Page coordinates: shift everything so the drawing starts one inch from the top left.
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  const take = (p: Point) => { x0 = Math.min(x0, p.x); y0 = Math.min(y0, p.y); x1 = Math.max(x1, p.x); y1 = Math.max(y1, p.y); };
  doc.atoms.forEach(take);
  doc.arrows.forEach((a) => { take({ x: a.x1, y: a.y1 }); take({ x: a.x2, y: a.y2 }); });
  doc.texts.forEach((t) => { const b = textBox(t); take(b); take({ x: b.x + b.w, y: b.y + b.h }); });
  if (x0 === Infinity) { x0 = y0 = x1 = y1 = 0; }
  const X = (x: number) => (x - x0) * PT + 72, Y = (y: number) => (y - y0) * PT + 72;
  const P = (x: number, y: number) => `${r2(X(x))} ${r2(Y(y))}`;

  const body: string[] = [];
  const frags = fragments(doc);
  const fragId = new Map<number, number>(); // atom id -> fragment id
  const brackets: string[] = [];
  const symbols: string[] = [];

  for (const comp of frags) {
    const fid = newId();
    comp.forEach((id) => fragId.set(id, fid));
    const inComp = new Set(comp);
    const out: string[] = [`<fragment id="${fid}">`];
    for (const id of comp) {
      const a = byId.get(id)!;
      const attrs = [`id="${a.id}"`, `p="${P(a.x, a.y)}"`];
      let label = '';
      let right = false;
      if (a.el === '*') {
        if (!ends.has(a.id)) {
          attrs.push('NodeType="GenericNickname"', 'GenericNickname="*"');
          label = '*';
        }
        // a repeat-unit end is a plain bond end crossing the bracket, as ChemDraw draws it
      } else {
        const z = Molecule.getAtomicNoFromLabel(a.el);
        if (z !== 6) attrs.push(`Element="${z}"`);
        if (a.charge) attrs.push(`Charge="${a.charge}"`);
        if (a.iso) attrs.push(`Isotope="${a.iso}"`);
        const li = labelInfo(doc, a, nbrs.get(a.id)!, an);
        if (a.radical) {
          attrs.push(`Radical="${a.radical[0].toUpperCase()}${a.radical.slice(1)}"`);
          // ChemDraw draws the dot as a symbol that represents the node's Radical attribute
          const dy = li.shown ? FONT * 0.62 : 6.5;
          symbols.push(`<graphic id="${newId()}" BoundingBox="${P(a.x, a.y - dy)} ${P(a.x, a.y - dy - 8)}" GraphicType="Symbol" SymbolType="${a.radical === 'doublet' ? 'Electron' : 'LonePair'}"><represent attribute="Radical" object="${a.id}"/></graphic>`);
        }
        if (li.shown) {
          attrs.push(`NumHydrogens="${li.hCount}"`);
          const h = li.hCount ? `H${li.hCount > 1 ? li.hCount : ''}` : '';
          right = !!h && li.hSide === 'left';
          const ch = a.charge ? `${Math.abs(a.charge) > 1 ? Math.abs(a.charge) : ''}${a.charge > 0 ? '+' : '-'}` : '';
          label = right ? `${h}${a.el}${ch}` : `${a.el}${h}${ch}`;
        }
      }
      if (!label) { out.push(`<n ${attrs.join(' ')}/>`); continue; }
      const half = textWidth(a.el === '*' ? '*' : a.el, LABEL_SIZE) / 2;
      const tp = right ? `${r2(X(a.x) + half)} ${r2(Y(a.y) + 3.5)}` : `${r2(X(a.x) - half)} ${r2(Y(a.y) + 3.5)}`;
      const just = right ? ' Justification="Right" LabelJustification="Right"' : ' LabelJustification="Left"';
      out.push(`<n ${attrs.join(' ')}><t p="${tp}"${just}><s font="3" size="${LABEL_SIZE}" face="96">${xmlEsc(label)}</s></t></n>`);
    }
    for (const b of doc.bonds) {
      if (!inComp.has(b.a)) continue;
      const attrs = [`id="${b.id}"`, `B="${b.a}"`, `E="${b.b}"`];
      if (b.order !== 1) attrs.push(`Order="${b.order}"`);
      if (b.order === 1 && b.stereo !== 'none') attrs.push(`Display="${b.stereo === 'wedge' ? 'WedgeBegin' : b.stereo === 'hash' ? 'WedgedHashBegin' : 'Wavy'}"`);
      if (b.order === 2 && b.stereo === 'wavy') attrs.push('Display="Wavy"');
      out.push(`<b ${attrs.join(' ')}/>`);
    }
    out.push('</fragment>');
    body.push(out.join(''));

  }

  // Repeat units: one bracket across each crossed bond, grouped as a structural repeat unit (SRU).
  // Automatic ones sit on the bonds to two * ends; Bracket-tool groups (copolymer blocks) inside a chain.
  const marks = bracketMarks(doc);
  for (const gid0 of new Set(marks.map((m) => m.group))) {
    const ms = marks.filter((m) => m.group === gid0);
    const group = (doc.brackets ?? []).find((g) => g.id === gid0);
    const label = ms.find((m) => m.sub)?.sub?.text ?? 'n';
    const inside = group ? group.atoms : (frags.find((c) => c.includes(gid0)) ?? []).filter((id) => !ends.has(id));
    const attach: string[] = [];
    for (const m of ms) {
      const gid = newId();
      const tag = m.sub
        ? `<objecttag TagType="String" Name="bracketusage" Value="${xmlEsc(label)}"><t p="${P(m.sub.x, m.sub.y)}"><s font="3" size="${LABEL_SIZE * 0.75}" face="2">${xmlEsc(label)}</s></t></objecttag>`
        : '';
      // the bracket's lips point to the left of its start-to-end direction, so run it along (uy, -ux)
      brackets.push(`<graphic id="${gid}" BoundingBox="${P(m.e1.x, m.e1.y)} ${P(m.e2.x, m.e2.y)}" GraphicType="Bracket" BracketType="Square" LipSize="60"${tag ? ' BracketUsage="SRU"' : ''}>${tag}</graphic>`);
      attach.push(`<bracketattachment id="${newId()}" GraphicID="${gid}"><crossingbond id="${newId()}" BondID="${m.bond}" InnerAtomID="${m.inner}"/></bracketattachment>`);
    }
    brackets.push(`<bracketedgroup id="${newId()}" BracketedObjectIDs="${inside.join(' ')}" BracketUsage="SRU" PolymerRepeatPattern="HeadToTail" SRULabel="${xmlEsc(label)}">${attach.join('')}</bracketedgroup>`);
  }

  for (const ar of doc.arrows) {
    body.push(`<graphic id="${ar.id}" BoundingBox="${P(ar.x2, ar.y2)} ${P(ar.x1, ar.y1)}" GraphicType="Line" ArrowType="${ARROW_TYPE[ar.kind]}" HeadSize="1000"/>`);
  }
  for (const t of doc.texts) {
    body.push(`<t id="${t.id}" p="${P(t.x, t.y)}" Justification="Left"><s font="3" size="${LABEL_SIZE}" face="96">${xmlEsc(t.text)}</s></t>`);
  }
  body.push(...brackets, ...symbols);
  body.push(schemeXml(doc, frags, fragId, newId));

  const W = Math.max(540, Math.ceil(X(x1) + 72)), H = Math.max(720, Math.ceil(Y(y1) + 72));
  return [
    '<?xml version="1.0" encoding="UTF-8" ?>',
    '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >',
    `<CDXML CreationProgram="MolCanvas" Name="MolCanvas drawing" BoundingBox="${r2(72)} ${r2(72)} ${r2(X(x1))} ${r2(Y(y1))}"`
      + ' InterpretChemically="yes" LabelFont="3" LabelSize="10" LabelFace="96" CaptionFont="3" CaptionSize="10"'
      + ' HashSpacing="2.50" MarginWidth="1.60" LineWidth="0.60" BoldWidth="2" BondLength="14.40" BondSpacing="18"'
      + ' ChainAngle="120" LabelJustification="Auto" CaptionJustification="Left" color="0" bgcolor="1">',
    '<colortable><color r="1" g="1" b="1"/><color r="0" g="0" b="0"/></colortable>',
    '<fonttable><font id="3" charset="iso-8859-1" name="Arial"/></fonttable>',
    `<page id="${newId()}" BoundingBox="0 0 ${W} ${H}" HeaderPosition="36" FooterPosition="36" HeightPages="1" WidthPages="1">`,
    ...body.filter(Boolean),
    '</page>',
    '</CDXML>',
    '',
  ].join('\n');
}

/**
 * Reaction steps for forward arrows, so ChemDraw (and other readers) know which structures are
 * reactants and products: the molecules in line with an arrow on its tail side are reactants, those
 * on its head side are products, up to the next arrow on the same line. Text over or under an arrow
 * is recorded as conditions.
 */
function schemeXml(doc: Doc, frags: number[][], fragId: Map<number, number>, newId: () => number): string {
  const steps: string[] = [];
  const byId = new Map(doc.atoms.map((a) => [a.id, a]));
  const centres = frags.map((c) => {
    const pts = c.map((id) => byId.get(id)!);
    return { id: fragId.get(c[0])!, x: pts.reduce((s, p) => s + p.x, 0) / pts.length, y: pts.reduce((s, p) => s + p.y, 0) / pts.length };
  });
  const forward = doc.arrows.filter((a) => a.kind === 'forward');
  for (const ar of forward) {
    const len = Math.hypot(ar.x2 - ar.x1, ar.y2 - ar.y1) || 1;
    const ux = (ar.x2 - ar.x1) / len, uy = (ar.y2 - ar.y1) / len;
    const along = (p: Point) => (p.x - ar.x1) * ux + (p.y - ar.y1) * uy;
    const across = (p: Point) => (p.x - ar.x1) * -uy + (p.y - ar.y1) * ux;
    const inLine = (p: Point) => Math.abs(across(p)) < BOND * 2.5;
    // other arrows on the same line limit how far reactants and products reach
    let lo = -Infinity, hi = Infinity;
    for (const o of doc.arrows) {
      if (o === ar) continue;
      const s = { x: o.x1, y: o.y1 }, e = { x: o.x2, y: o.y2 };
      if (!inLine(s) || !inLine(e)) continue;
      const a0 = Math.min(along(s), along(e)), a1 = Math.max(along(s), along(e));
      if (a1 <= 0) lo = Math.max(lo, a1);
      if (a0 >= len) hi = Math.min(hi, a0);
    }
    const reactants = centres.filter((c) => inLine(c) && along(c) < 0 && along(c) > lo).map((c) => c.id);
    const products = centres.filter((c) => inLine(c) && along(c) > len && along(c) < hi).map((c) => c.id);
    if (!reactants.length || !products.length) continue;
    const near = (t: TextItem) => {
      const b = textBox(t), c = { x: b.x + b.w / 2, y: b.y + b.h / 2 };
      return along(c) > 0 && along(c) < len && Math.abs(across(c)) < BOND * 1.5 ? across(c) : null;
    };
    const above = doc.texts.filter((t) => (near(t) ?? 0) < 0).map((t) => t.id);
    const below = doc.texts.filter((t) => (near(t) ?? 0) > 0).map((t) => t.id);
    steps.push(`<step id="${newId()}" ReactionStepReactants="${reactants.join(' ')}" ReactionStepProducts="${products.join(' ')}" ReactionStepArrows="${ar.id}"`
      + `${above.length ? ` ReactionStepObjectsAboveArrow="${above.join(' ')}"` : ''}${below.length ? ` ReactionStepObjectsBelowArrow="${below.join(' ')}"` : ''}/>`);
  }
  return steps.length ? `<scheme id="${newId()}">${steps.join('')}</scheme>` : '';
}

