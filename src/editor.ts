// Interactive editor: tools, hit testing, selection, undo/redo, pan/zoom.
import { analyze, Analysis, cleanCoordinates, fromMolfile, fromSmiles, toSmiles } from './chem';
import { fromCdxml } from './cdxml';
import {
  addAtom, addBond, Arrow, ArrowKind, atomById, BOND, BondStereo, bondsOf, cloneDoc, Doc, emptyDoc, fragments,
  mergeInto, neighbors, Point, Radical, bounds, pruneLoneCarbons, removeObjects, segDist, subDoc, dist, centroid, TextItem,
  bracketProblem, repeatUnitEnds,
} from './model';
import {
  chain, dirVec, findAtomNear, freeDirection, Placement, placeTemplate, ringAtPoint, ringOnAtom, ringOnBond,
  orientPolymer, rotatePoint, singleBond, snapAngle,
} from './geometry';
import { bracketMarks, renderDoc, RenderOpts, textBox, Theme } from './render';
import type { Preset } from './presets';
import { Explanation, explainAtom, explainBond, explainProblem, valenceProblems } from './valence';

export type Tool =
  | { kind: 'select' }
  | { kind: 'erase' }
  | { kind: 'bond'; order: 1 | 2 | 3; stereo: BondStereo }
  | { kind: 'chain' }
  | { kind: 'ring'; size: number; aromatic: boolean }
  | { kind: 'atom'; el: string }
  | { kind: 'charge'; delta: 1 | -1 }
  | { kind: 'radical' }
  | { kind: 'arrow'; arrow: ArrowKind }
  | { kind: 'text' }
  | { kind: 'bracket' }
  | { kind: 'template'; preset: Preset; frag: Doc };

type Hit =
  | { type: 'atom'; id: number }
  | { type: 'bond'; id: number }
  | { type: 'arrow'; id: number }
  | { type: 'text'; id: number }
  | { type: 'bracket'; id: number }
  | null;

interface Drag {
  start: Point; // world
  startScreen: Point;
  hit: Hit;
  moved: boolean;
  mode: 'move' | 'rubber' | 'draw' | 'pan' | 'erase' | 'arrow' | 'none';
  original?: Doc;
  moveIds?: Set<number>;
  shift: boolean;
}

const HOTKEY_ELEMENTS: Record<string, string> = {
  c: 'C', n: 'N', o: 'O', s: 'S', p: 'P', f: 'F', h: 'H', i: 'I', l: 'Cl', b: 'Br', k: 'K', '*': '*',
  C: 'C', N: 'N', O: 'O', S: 'S', P: 'P', F: 'F', H: 'H', I: 'I', L: 'Cl', B: 'B',
};

export class Editor {
  doc: Doc = emptyDoc();
  sel = new Set<number>();
  tool: Tool = { kind: 'bond', order: 1, stereo: 'none' };
  hover: Hit = null;
  zoom = 1;
  tx = 0;
  ty = 0;
  theme!: Theme;
  colorAtoms = true;
  showLonePairs = false;
  showErrors = true;

  onChange: () => void = () => {};
  onStatus: (msg: string) => void = () => {};
  onToolChange: (t: Tool) => void = () => {};
  /** Explanation of whatever is under the pointer (null: nothing, show the tool hint). */
  onHover: (info: Explanation | null) => void = () => {};
  /** A change was refused because it would break a chemistry rule. */
  onWarn: (msg: string) => void = () => {};
  /** Set when the previewed action would break a rule; shown in red and refused on click. */
  private ghostProblem: string | null = null;

  private undoStack: string[] = [];
  private redoStack: string[] = [];
  private analysisCache: { version: number; an: Analysis } | null = null;
  private version = 0;
  private drag: Drag | null = null;
  private ghost: Doc | null = null;
  private pointer: Point = { x: 0, y: 0 };
  private world: SVGGElement;
  private spaceDown = false;
  private input: HTMLInputElement | null = null;
  private clipboard: { json: string; smiles: string } | null = null;

  constructor(private svg: SVGSVGElement, private host: HTMLElement, theme: Theme) {
    this.theme = theme;
    this.world = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    svg.appendChild(this.world);
    this.bindEvents();
    requestAnimationFrame(() => this.centerView());
  }

  // ---------- document state ----------

  private analysis(): Analysis {
    if (!this.analysisCache || this.analysisCache.version !== this.version) {
      this.analysisCache = { version: this.version, an: analyze(this.doc) };
    }
    return this.analysisCache.an;
  }

  /**
   * Apply a change with undo support. Changes that would push an atom past its valence
   * (C with 5 bonds, a polymer end with 2…) are refused and explained via onWarn.
   */
  commit(mutate: (d: Doc) => void): boolean {
    const before = JSON.stringify(this.doc);
    const prior = new Map(valenceProblems(this.doc).map((p) => [p.atom.id, p.used]));
    mutate(this.doc);
    const after = JSON.stringify(this.doc);
    if (before === after) return false;
    const fresh = valenceProblems(this.doc).filter((p) => p.used > (prior.get(p.atom.id) ?? -1));
    if (fresh.length) {
      this.doc = JSON.parse(before);
      this.version++;
      this.render();
      this.onWarn(explainProblem(fresh[0]));
      return false;
    }
    this.undoStack.push(before);
    if (this.undoStack.length > 300) this.undoStack.shift();
    this.redoStack = [];
    this.changed();
    return true;
  }

  /**
   * Radical tool: none → • (doublet) → •• (triplet) → •• (singlet) → none. States the atom has no
   * electrons left for are skipped (CH3• can become :CH2 only if it loses an H, so it goes back to none).
   */
  cycleRadical(id: number) {
    const a = atomById(this.doc, id);
    if (!a) return;
    if (a.el === '*') { this.onWarn('A polymer end (*) marks a bond position, not an atom, so it cannot carry a radical.'); return; }
    const order: (Radical | undefined)[] = [undefined, 'doublet', 'triplet', 'singlet'];
    const set = (r: Radical | undefined) => (d: Doc) => { const x = atomById(d, id)!; if (r) x.radical = r; else delete x.radical; };
    const start = order.indexOf(a.radical);
    // from no radical, try the doublet directly so a refusal is explained in the info bar
    if (start === 0) { this.commit(set('doublet')); return; }
    const prior = new Map(valenceProblems(this.doc).map((p) => [p.atom.id, p.used]));
    for (let k = start + 1; k < order.length; k++) {
      const trial = cloneDoc(this.doc);
      set(order[k])(trial);
      if (!valenceProblems(trial).some((p) => p.used > (prior.get(p.atom.id) ?? -1))) { this.commit(set(order[k])); return; }
    }
    this.commit(set(undefined));
  }

  /** Why a placement would be refused, or null if it is chemically fine. */
  private placementProblem(pl: Placement): string | null {
    const prior = new Map(valenceProblems(this.doc).map((p) => [p.atom.id, p.used]));
    const d = cloneDoc(this.doc);
    this.applyPlacement(d, pl);
    const fresh = valenceProblems(d).filter((p) => p.used > (prior.get(p.atom.id) ?? -1));
    return fresh.length ? explainProblem(fresh[0]) : null;
  }

  private setGhost(pl: Placement | null) {
    this.ghost = pl ? this.ghostOf(pl) : null;
    this.ghostProblem = pl ? this.placementProblem(pl) : null;
  }

  /** Plain-language description of an atom or bond for the info bar. */
  describe(hit: Hit): Explanation | null {
    if (!hit) return null;
    const an = this.analysis();
    if (hit.type === 'atom') {
      const a = atomById(this.doc, hit.id)!;
      const inAromatic = this.doc.bonds.some((b) => (b.a === a.id || b.b === a.id) && an.ringCenter.has(b.id));
      return explainAtom(this.doc, a, an.hCount.get(a.id) ?? 0, inAromatic);
    }
    if (hit.type === 'bond') {
      const b = this.doc.bonds.find((x) => x.id === hit.id)!;
      return explainBond(this.doc, b, an.ringCenter.has(b.id));
    }
    if (hit.type === 'arrow') {
      const ar = this.doc.arrows.find((x) => x.id === hit.id)!;
      const what = { forward: 'Reaction arrow: reactants on the left give products on the right.', equilibrium: 'Equilibrium arrows: the reaction runs in both directions.', retro: 'Retrosynthesis arrow: the target (left) can be made from the precursors (right).', resonance: 'Resonance arrow: both structures are contributors to one real molecule.' }[ar.kind];
      return { title: 'Arrow', detail: `${what} Drag with the Select tool to move it.` };
    }
    if (hit.type === 'bracket') {
      const g = (this.doc.brackets ?? []).find((x) => x.id === hit.id)!;
      return { title: `Repeat unit ×${g.sub}`, detail: `The bracketed block repeats ${g.sub} times along the chain. Each bracketed block keeps its own count, so a block copolymer reads –[A]ₙ–[B]ₘ–. Click it with the Bracket or Text tool to change the count; erase it to remove the brackets.` };
    }
    const t = this.doc.texts.find((x) => x.id === hit.id)!;
    return { title: 'Text', detail: `“${t.text}”. Double-click to edit. Digits after letters are shown as subscripts.` };
  }

  private emitHover() {
    if (this.ghostProblem) this.onHover({ title: 'Not allowed', detail: this.ghostProblem, warn: true });
    else this.onHover(this.describe(this.hover));
  }

  private changed() {
    this.version++;
    const alive = new Set<number>([...this.doc.atoms, ...this.doc.bonds, ...this.doc.arrows, ...this.doc.texts].map((o) => o.id));
    for (const id of [...this.sel]) if (!alive.has(id)) this.sel.delete(id);
    this.render();
    this.onChange();
  }

  undo() {
    const prev = this.undoStack.pop();
    if (!prev) return;
    this.redoStack.push(JSON.stringify(this.doc));
    this.doc = JSON.parse(prev);
    this.changed();
  }

  redo() {
    const next = this.redoStack.pop();
    if (!next) return;
    this.undoStack.push(JSON.stringify(this.doc));
    this.doc = JSON.parse(next);
    this.changed();
  }

  canUndo() { return this.undoStack.length > 0; }
  canRedo() { return this.redoStack.length > 0; }

  /** Replace the drawing. False if it was refused because it breaks a valence rule. */
  load(doc: Doc, keepHistory = true): boolean {
    let ok = true;
    if (keepHistory) ok = this.commit((d) => Object.assign(d, cloneDoc(doc)));
    else { this.doc = cloneDoc(doc); this.undoStack = []; this.redoStack = []; this.changed(); }
    this.sel.clear();
    this.fitView();
    return ok;
  }

  clear() {
    this.commit((d) => Object.assign(d, emptyDoc(), { nextId: d.nextId }));
  }

  setTool(t: Tool) {
    this.tool = t;
    this.pendingBracketBond = null;
    this.ghost = null;
    this.ghostProblem = null;
    this.closeInput();
    this.svg.dataset.tool = t.kind;
    this.onToolChange(t);
    this.onStatus(this.toolHint());
    this.render();
  }

  toolHint(): string {
    const t = this.tool;
    switch (t.kind) {
      case 'select': return 'Select: drag to box-select, drag a selection to move it, double-click to select a whole molecule.';
      case 'erase': return 'Eraser: click or drag over atoms, bonds, arrows and text to delete them.';
      case 'bond': return 'Bond: click an atom to grow a bond, drag to set direction, click a bond to change it. Hover + letter key sets an element.';
      case 'chain': return 'Chain: drag from an atom or empty space to draw a zig-zag carbon chain.';
      case 'ring': return 'Ring: click empty space to place, click a bond to fuse, click an atom for a spiro/substituent ring.';
      case 'atom': return `Atom ${t.el}: click an atom to change it, click empty space to add one, drag from an atom to add a bonded ${t.el}.`;
      case 'charge': return `Charge: click an atom to ${t.delta > 0 ? 'add a positive' : 'add a negative'} charge (or cancel the opposite one).`;
      case 'radical': return 'Radical: click an atom to cycle • (radical) → •• triplet → •• singlet → none. Hover an atom and press . for the same.';
      case 'arrow': return 'Arrow: drag to draw a reaction arrow (snaps to 15°), click to drop a standard arrow.';
      case 'bracket': return 'Bracket: click a bond on one side of a repeat unit, then a bond on the other side (or drag a box around it). Click a bracket to set its count; each block of a copolymer has its own.';
      case 'text': return 'Text: click to add a label (digits after letters become subscripts, e.g. CH2Cl2). Click a polymer end to edit n.';
      case 'template': return `${t.preset.name}: click to place${t.preset.kind === 'group' ? ', or click an atom to attach' : ''}. Esc to stop.`;
    }
  }

  // ---------- coordinates ----------

  toWorld(clientX: number, clientY: number): Point {
    const r = this.svg.getBoundingClientRect();
    return { x: (clientX - r.left - this.tx) / this.zoom, y: (clientY - r.top - this.ty) / this.zoom };
  }

  viewCenter(): Point {
    const r = this.svg.getBoundingClientRect();
    return { x: (r.width / 2 - this.tx) / this.zoom, y: (r.height / 2 - this.ty) / this.zoom };
  }

  centerView() {
    const r = this.svg.getBoundingClientRect();
    this.tx = r.width / 2;
    this.ty = r.height / 2;
    this.render();
  }

  fitView() {
    const r = this.svg.getBoundingClientRect();
    if (!this.doc.atoms.length && !this.doc.arrows.length && !this.doc.texts.length) {
      this.zoom = 1;
      this.centerView();
      return;
    }
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    const pts: Point[] = [...this.doc.atoms, ...this.doc.arrows.flatMap((a) => [{ x: a.x1, y: a.y1 }, { x: a.x2, y: a.y2 }]), ...this.doc.texts];
    for (const p of pts) { x0 = Math.min(x0, p.x); y0 = Math.min(y0, p.y); x1 = Math.max(x1, p.x); y1 = Math.max(y1, p.y); }
    const w = x1 - x0 + BOND * 3, h = y1 - y0 + BOND * 3;
    this.zoom = Math.max(0.2, Math.min(1.6, Math.min(r.width / w, r.height / h)));
    this.tx = r.width / 2 - ((x0 + x1) / 2) * this.zoom;
    this.ty = r.height / 2 - ((y0 + y1) / 2) * this.zoom;
    this.render();
  }

  zoomBy(factor: number, around?: Point) {
    const r = this.svg.getBoundingClientRect();
    const sx = around ? around.x : r.width / 2, sy = around ? around.y : r.height / 2;
    const nz = Math.max(0.15, Math.min(6, this.zoom * factor));
    this.tx = sx - ((sx - this.tx) * nz) / this.zoom;
    this.ty = sy - ((sy - this.ty) * nz) / this.zoom;
    this.zoom = nz;
    this.render();
  }

  // ---------- rendering ----------

  render() {
    const opts: RenderOpts = {
      theme: this.theme,
      colorAtoms: this.colorAtoms,
      selection: this.sel,
      hover: this.hover ? this.hover.id : null,
      analysis: this.analysis(),
      showErrors: this.showErrors,
      showLonePairs: this.showLonePairs,
    };
    let html = renderDoc(this.doc, opts);
    if (this.ghost) {
      const theme = this.ghostProblem ? { ...this.theme, fg: this.theme.error } : this.theme;
      html += `<g opacity="${this.ghostProblem ? 0.6 : 0.45}" pointer-events="none">${renderDoc(this.ghost, { theme, colorAtoms: this.colorAtoms && !this.ghostProblem })}</g>`;
    }
    if (this.drag?.mode === 'rubber' && this.drag.moved) {
      const p = this.pointer, s = this.drag.start;
      html += `<rect x="${Math.min(p.x, s.x)}" y="${Math.min(p.y, s.y)}" width="${Math.abs(p.x - s.x)}" height="${Math.abs(p.y - s.y)}" fill="${this.theme.select}" fill-opacity="0.08" stroke="${this.theme.select}" stroke-width="${1 / this.zoom}" stroke-dasharray="${4 / this.zoom} ${3 / this.zoom}"/>`;
    }
    const pb = this.tool.kind === 'bracket' && this.pendingBracketBond !== null ? this.doc.bonds.find((b) => b.id === this.pendingBracketBond) : undefined;
    if (pb) {
      const A = atomById(this.doc, pb.a)!, B = atomById(this.doc, pb.b)!;
      html += `<line x1="${A.x}" y1="${A.y}" x2="${B.x}" y2="${B.y}" stroke="${this.theme.select}" stroke-width="10" stroke-opacity="0.45" stroke-linecap="round" pointer-events="none"/>`;
    }
    this.world.innerHTML = html;
    this.world.setAttribute('transform', `translate(${this.tx},${this.ty}) scale(${this.zoom})`);
    const g = 20 * this.zoom;
    this.host.style.backgroundSize = `${g}px ${g}px`;
    this.host.style.backgroundPosition = `${this.tx}px ${this.ty}px`;
  }

  // ---------- hit testing ----------

  hitTest(p: Point): Hit {
    const tol = Math.max(6, 9 / this.zoom);
    let best: Hit = null, bd = tol;
    for (const a of this.doc.atoms) {
      const d = dist(a, p);
      if (d < bd) { bd = d; best = { type: 'atom', id: a.id }; }
    }
    if (best) return best;
    const bh = this.bracketHit(p);
    if (bh) return bh;
    bd = Math.max(4, 6 / this.zoom);
    for (const b of this.doc.bonds) {
      const A = atomById(this.doc, b.a)!, B = atomById(this.doc, b.b)!;
      const d = segDist(p, A, B);
      if (d < bd) { bd = d; best = { type: 'bond', id: b.id }; }
    }
    if (best) return best;
    for (const t of this.doc.texts) {
      const bx = textBox(t);
      if (p.x >= bx.x && p.x <= bx.x + bx.w && p.y >= bx.y && p.y <= bx.y + bx.h) return { type: 'text', id: t.id };
    }
    for (const ar of this.doc.arrows) {
      if (segDist(p, { x: ar.x1, y: ar.y1 }, { x: ar.x2, y: ar.y2 }) < 7) return { type: 'arrow', id: ar.id };
    }
    return null;
  }

  /** A bracket line or its count label under p (in-chain bracket groups; * end counts via the Text tool). */
  private bracketHit(p: Point): Hit {
    const tol = Math.max(4, 6 / this.zoom);
    const groups = new Set((this.doc.brackets ?? []).map((g) => g.id));
    const k = this.tool.kind;
    for (const m of bracketMarks(this.doc)) {
      const own = groups.has(m.group);
      if (!own && k !== 'text' && k !== 'bracket') continue;
      const onSub = m.sub && p.x >= m.sub.x - 3 && p.x <= m.sub.x + 8 * m.sub.text.length + 3 && p.y >= m.sub.y - 12 && p.y <= m.sub.y + 3;
      // with the Bracket tool, a click on the bond itself picks the bond (to start or end a unit)
      const onBond = k === 'bracket' && (() => { const b = this.doc.bonds.find((x) => x.id === m.bond)!; return segDist(p, atomById(this.doc, b.a)!, atomById(this.doc, b.b)!) < tol; })();
      if (onSub || (!onBond && segDist(p, m.e1, m.e2) < tol)) return own ? { type: 'bracket', id: m.group } : { type: 'atom', id: m.group };
    }
    return null;
  }

  // ---------- events ----------

  private bindEvents() {
    const svg = this.svg;
    svg.addEventListener('pointerdown', (e) => this.pointerDown(e));
    svg.addEventListener('pointermove', (e) => this.pointerMove(e));
    svg.addEventListener('pointerup', (e) => this.pointerUp(e));
    svg.addEventListener('pointerleave', () => { if (!this.drag) { this.hover = null; this.ghost = null; this.ghostProblem = null; this.render(); this.emitHover(); } });
    svg.addEventListener('dblclick', (e) => this.doubleClick(e));
    svg.addEventListener('wheel', (e) => {
      e.preventDefault();
      if (e.ctrlKey || e.metaKey) {
        const r = svg.getBoundingClientRect();
        this.zoomBy(Math.exp(-e.deltaY * 0.01), { x: e.clientX - r.left, y: e.clientY - r.top });
      } else {
        this.tx -= e.deltaX;
        this.ty -= e.deltaY;
        this.render();
      }
    }, { passive: false });
    window.addEventListener('keydown', (e) => this.keyDown(e));
    window.addEventListener('keyup', (e) => { if (e.code === 'Space') { this.spaceDown = false; svg.style.cursor = ''; } });
    document.addEventListener('copy', (e) => this.copy(e, false));
    document.addEventListener('cut', (e) => this.copy(e, true));
    document.addEventListener('paste', (e) => this.paste(e));
    new ResizeObserver(() => this.render()).observe(svg);
  }

  private pointerDown(e: PointerEvent) {
    if (this.input) { this.commitInput(); return; }
    const p = this.toWorld(e.clientX, e.clientY);
    this.svg.setPointerCapture(e.pointerId);
    const hit = this.hitTest(p);
    const drag: Drag = { start: p, startScreen: { x: e.clientX, y: e.clientY }, hit, moved: false, mode: 'none', shift: e.shiftKey };
    this.drag = drag;
    if (e.button === 1 || this.spaceDown) { drag.mode = 'pan'; return; }
    const t = this.tool;
    if (t.kind === 'select') {
      if (hit) {
        if (!this.sel.has(hit.id)) {
          if (!e.shiftKey) this.sel.clear();
          this.selectHit(hit);
        } else if (e.shiftKey) {
          this.sel.delete(hit.id);
        }
        drag.mode = 'move';
        drag.original = cloneDoc(this.doc);
        drag.moveIds = this.movableAtoms();
      } else {
        if (!e.shiftKey) this.sel.clear();
        drag.mode = 'rubber';
      }
      this.render();
    } else if (t.kind === 'erase') {
      drag.mode = 'erase';
      drag.original = cloneDoc(this.doc);
      this.eraseAt(p, true);
    } else if (t.kind === 'bracket') {
      drag.mode = hit ? 'none' : 'rubber';
    } else if (t.kind === 'arrow') {
      drag.mode = 'arrow';
    } else if (t.kind === 'bond' || t.kind === 'chain' || t.kind === 'atom') {
      drag.mode = 'draw';
    } else {
      drag.mode = 'none';
    }
  }

  private pointerMove(e: PointerEvent) {
    const p = this.toWorld(e.clientX, e.clientY);
    this.pointer = p;
    const drag = this.drag;
    if (drag && !drag.moved && Math.hypot(e.clientX - drag.startScreen.x, e.clientY - drag.startScreen.y) > 4) drag.moved = true;
    if (drag?.mode === 'pan') {
      this.tx += e.movementX;
      this.ty += e.movementY;
      this.render();
      return;
    }
    if (drag && drag.moved) {
      if (drag.mode === 'move' && drag.original && drag.moveIds) {
        let dx = p.x - drag.start.x, dy = p.y - drag.start.y;
        const objs = this.sel;
        this.doc = cloneDoc(drag.original);
        for (const a of this.doc.atoms) if (drag.moveIds.has(a.id)) { a.x += dx; a.y += dy; }
        for (const ar of this.doc.arrows) if (objs.has(ar.id)) { ar.x1 += dx; ar.y1 += dy; ar.x2 += dx; ar.y2 += dy; }
        for (const t of this.doc.texts) if (objs.has(t.id)) { t.x += dx; t.y += dy; }
        this.version++;
        this.render();
        return;
      }
      if (drag.mode === 'rubber') { this.render(); return; }
      if (drag.mode === 'erase') { this.eraseAt(p, false); return; }
      if (drag.mode === 'arrow') {
        const ar = this.arrowFromDrag(drag.start, p);
        this.ghost = { ...emptyDoc(), arrows: [ar] };
        this.render();
        return;
      }
      if (drag.mode === 'draw') {
        const pl = this.drawPlacement(drag, p);
        this.setGhost(pl);
        if (this.ghostProblem) this.emitHover();
        else if (this.tool.kind === 'chain' && pl && 'count' in pl) this.onHover({ title: 'Chain', detail: `${(pl as Placement & { count: number }).count} bonds. Release to add the chain.` });
        this.render();
        return;
      }
    }
    if (!drag) {
      const hit = this.hitTest(p);
      const changed = JSON.stringify(hit) !== JSON.stringify(this.hover);
      this.hover = hit;
      const t = this.tool.kind;
      const previewOnAtom = (t === 'bond' || t === 'atom' || t === 'chain') && hit?.type === 'atom';
      if (t === 'ring' || t === 'template' || previewOnAtom) {
        this.setGhost(this.clickPlacement(p, hit));
        this.render();
        if (changed || this.ghostProblem) this.emitHover();
      } else {
        if (this.ghost) { this.ghost = null; this.ghostProblem = null; this.render(); }
        else if (changed) this.render();
        if (changed) this.emitHover();
      }
    }
  }

  private pointerUp(e: PointerEvent) {
    const drag = this.drag;
    this.drag = null;
    if (!drag) return;
    const p = this.toWorld(e.clientX, e.clientY);
    this.ghost = null;
    this.ghostProblem = null;
    const t = this.tool;
    if (drag.mode === 'pan') return;
    if (drag.mode === 'move' && drag.original) {
      if (drag.moved) {
        const after = this.doc;
        this.doc = drag.original;
        this.commit((d) => {
          Object.assign(d, after);
          this.mergeDroppedAtom(d, drag);
        });
      }
      this.render();
      return;
    }
    if (drag.mode === 'rubber' && t.kind === 'bracket') {
      if (drag.moved) this.bracketBox(drag.start, p);
      else this.click(p, drag.hit, e);
      return;
    }
    if (drag.mode === 'rubber') {
      if (drag.moved) this.rubberSelect(drag.start, p);
      this.render();
      this.onChange();
      return;
    }
    if (drag.mode === 'erase' && drag.original) {
      const after = this.doc;
      this.doc = drag.original;
      this.commit((d) => Object.assign(d, after));
      return;
    }
    if (drag.mode === 'arrow') {
      const ar = drag.moved ? this.arrowFromDrag(drag.start, p) : null;
      if (!drag.moved && drag.hit?.type === 'arrow' && t.kind === 'arrow') {
        const id = drag.hit.id;
        this.commit((d) => { const a = d.arrows.find((x) => x.id === id); if (a) a.kind = t.arrow; });
        return;
      }
      this.commit((d) => {
        const a = ar ?? { id: 0, x1: p.x, y1: p.y, x2: p.x + BOND * 2.5, y2: p.y, kind: t.kind === 'arrow' ? t.arrow : 'forward' };
        a.id = d.nextId++;
        d.arrows.push(a);
      });
      return;
    }
    if (drag.mode === 'draw' && drag.moved) {
      const pl = this.drawPlacement(drag, p);
      if (pl) this.commit((d) => this.applyPlacement(d, pl));
      return;
    }
    // Plain click
    this.click(p, drag.hit, e);
  }

  private doubleClick(e: MouseEvent) {
    const p = this.toWorld(e.clientX, e.clientY);
    const hit = this.hitTest(p);
    if (this.tool.kind === 'select' && hit && (hit.type === 'atom' || hit.type === 'bond')) {
      const atomId = hit.type === 'atom' ? hit.id : this.doc.bonds.find((b) => b.id === hit.id)!.a;
      const comp = fragments(this.doc).find((c) => c.includes(atomId)) ?? [];
      for (const id of comp) this.sel.add(id);
      for (const b of this.doc.bonds) if (this.sel.has(b.a) && this.sel.has(b.b)) this.sel.add(b.id);
      this.render();
      this.onChange();
    } else if (hit?.type === 'text') {
      this.editText(hit.id);
    }
  }

  private click(p: Point, hit: Hit, e: PointerEvent) {
    const t = this.tool;
    switch (t.kind) {
      case 'bond': {
        if (hit?.type === 'bond') {
          const id = hit.id;
          this.commit((d) => {
            const b = d.bonds.find((x) => x.id === id)!;
            if (t.stereo !== 'none') {
              if (b.stereo === t.stereo) { const tmp = b.a; b.a = b.b; b.b = tmp; } // flip direction
              b.order = 1;
              b.stereo = t.stereo;
            } else if (t.order === 1) {
              b.stereo = 'none';
              b.order = b.order === 3 ? 1 : ((b.order + 1) as 1 | 2 | 3);
            } else {
              b.stereo = 'none';
              b.order = b.order === t.order ? 1 : t.order;
            }
          });
          return;
        }
        const pl = this.clickPlacement(p, hit);
        if (pl) this.commit((d) => this.applyPlacement(d, pl));
        return;
      }
      case 'chain':
      case 'ring':
      case 'template': {
        const pl = this.clickPlacement(p, hit);
        if (pl) this.commit((d) => this.applyPlacement(d, pl));
        return;
      }
      case 'atom': {
        if (hit?.type === 'atom') {
          const id = hit.id;
          this.commit((d) => { const a = atomById(d, id)!; a.el = t.el; delete a.sub; if (t.el === '*') delete a.radical; });
        } else if (!hit) {
          this.commit((d) => addAtom(d, t.el, p.x, p.y));
        }
        return;
      }
      case 'charge': {
        if (hit?.type === 'atom') {
          const id = hit.id;
          this.commit((d) => { const a = atomById(d, id)!; a.charge = Math.max(-4, Math.min(4, a.charge + t.delta)); });
        }
        return;
      }
      case 'radical': {
        if (hit?.type === 'atom') this.cycleRadical(hit.id);
        return;
      }
      case 'bracket': {
        if (hit?.type === 'bracket') this.editBracket(hit.id);
        else if (hit?.type === 'bond') this.bracketBondClick(hit.id);
        else if (hit?.type === 'atom' && atomById(this.doc, hit.id)?.el === '*') this.editAtom(hit.id);
        else this.onWarn('Click a bond on one side of the repeat unit, then a bond on the other side. You can also drag a box around the unit.');
        return;
      }
      case 'text': {
        if (hit?.type === 'bracket') this.editBracket(hit.id);
        else if (hit?.type === 'text') this.editText(hit.id);
        else if (hit?.type === 'atom') this.editAtom(hit.id);
        else if (!hit) this.newText(p);
        return;
      }
      case 'select': {
        if (!hit && !e.shiftKey) { this.sel.clear(); this.render(); this.onChange(); }
        return;
      }
    }
  }

  // ---------- placements ----------

  /** What a click at p (on `hit`) would add for the current tool. */
  private clickPlacement(p: Point, hit: Hit): Placement | null {
    const t = this.tool;
    if (t.kind === 'ring') {
      if (hit?.type === 'bond') return ringOnBond(this.doc, hit.id, t.size, t.aromatic, this.analysis().ringCenter.get(hit.id), p);
      if (hit?.type === 'atom') return ringOnAtom(this.doc, hit.id, t.size, t.aromatic);
      if (!hit) return ringAtPoint(this.doc, p, t.size, t.aromatic);
      return null;
    }
    if (t.kind === 'template') {
      if (hit?.type === 'atom' && t.preset.kind === 'group') return placeTemplate(this.doc, t.frag, 'group', p, hit.id);
      return placeTemplate(this.doc, t.frag, t.preset.kind, p);
    }
    if (t.kind === 'bond' || t.kind === 'chain' || t.kind === 'atom') {
      const order = t.kind === 'bond' ? t.order : 1;
      const stereo = t.kind === 'bond' ? t.stereo : 'none';
      const el = t.kind === 'atom' ? t.el : 'C';
      if (hit?.type === 'atom') {
        const A = atomById(this.doc, hit.id)!;
        const v = dirVec(freeDirection(this.doc, hit.id));
        return singleBond(this.doc, A, { x: A.x + v.x * BOND, y: A.y + v.y * BOND }, hit.id, el, order, stereo);
      }
      if (!hit && t.kind !== 'atom') {
        const v = dirVec(-Math.PI / 6);
        const s = { x: p.x - (v.x * BOND) / 2, y: p.y - (v.y * BOND) / 2 };
        return singleBond(this.doc, s, { x: s.x + v.x * BOND, y: s.y + v.y * BOND }, undefined, el, order, stereo);
      }
    }
    return null;
  }

  /** What dragging from drag.start to p would add. */
  private drawPlacement(drag: Drag, p: Point): Placement | null {
    const t = this.tool;
    const startAtom = drag.hit?.type === 'atom' ? drag.hit.id : undefined;
    const start = startAtom !== undefined ? atomById(this.doc, startAtom)! : drag.start;
    if (t.kind === 'chain') return chain(this.doc, start, p, startAtom);
    if (t.kind === 'bond' || t.kind === 'atom') {
      if (t.kind === 'atom' && startAtom === undefined) return null;
      const target = findAtomNear(this.doc, p, BOND * 0.3, new Set(startAtom !== undefined ? [startAtom] : []));
      let end: Point;
      if (target) end = target;
      else {
        const ang = snapAngle(Math.atan2(p.y - start.y, p.x - start.x), drag.shift ? 1 : 30);
        end = { x: start.x + Math.cos(ang) * BOND, y: start.y + Math.sin(ang) * BOND };
      }
      const order = t.kind === 'bond' ? t.order : 1;
      const stereo = t.kind === 'bond' ? t.stereo : 'none';
      return singleBond(this.doc, start, end, startAtom, t.kind === 'atom' ? t.el : 'C', order, stereo);
    }
    return null;
  }

  private applyPlacement(d: Doc, pl: Placement) {
    // Fragment ids are local; remap them to fresh ids in d.
    const newIds = mergeInto(d, pl.frag, 0, 0, pl.map);
    // Newly drawn atoms with a non-carbon element (atom tool) keep their element through mergeInto.
    void newIds;
  }

  private ghostOf(pl: Placement): Doc {
    // Show the fragment, with mapped atoms drawn at the existing atoms' positions.
    const g = cloneDoc(pl.frag);
    for (const a of g.atoms) {
      const ex = pl.map.get(a.id);
      if (ex !== undefined) {
        const A = atomById(this.doc, ex);
        if (A) { a.x = A.x; a.y = A.y; a.el = A.el === '*' ? 'C' : A.el; }
      }
    }
    // Hide mapped atoms' labels by making them carbon unless heteroatom (already drawn underneath).
    return g;
  }

  private arrowFromDrag(s: Point, p: Point): Arrow {
    const ang = snapAngle(Math.atan2(p.y - s.y, p.x - s.x), 15);
    const len = Math.max(BOND, dist(s, p));
    const t = this.tool;
    return { id: 0, x1: s.x, y1: s.y, x2: s.x + Math.cos(ang) * len, y2: s.y + Math.sin(ang) * len, kind: t.kind === 'arrow' ? t.arrow : 'forward' };
  }

  // ---------- selection ----------

  private selectHit(hit: NonNullable<Hit>) {
    this.sel.add(hit.id);
    if (hit.type === 'bond') {
      const b = this.doc.bonds.find((x) => x.id === hit.id)!;
      this.sel.add(b.a);
      this.sel.add(b.b);
    }
  }

  private movableAtoms(): Set<number> {
    const ids = new Set<number>();
    for (const a of this.doc.atoms) if (this.sel.has(a.id)) ids.add(a.id);
    for (const b of this.doc.bonds) if (this.sel.has(b.id)) { ids.add(b.a); ids.add(b.b); }
    return ids;
  }

  private rubberSelect(a: Point, b: Point) {
    const x0 = Math.min(a.x, b.x), x1 = Math.max(a.x, b.x), y0 = Math.min(a.y, b.y), y1 = Math.max(a.y, b.y);
    const inside = (p: Point) => p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1;
    for (const at of this.doc.atoms) if (inside(at)) this.sel.add(at.id);
    for (const bd of this.doc.bonds) if (this.sel.has(bd.a) && this.sel.has(bd.b)) this.sel.add(bd.id);
    for (const ar of this.doc.arrows) if (inside({ x: ar.x1, y: ar.y1 }) && inside({ x: ar.x2, y: ar.y2 })) this.sel.add(ar.id);
    for (const t of this.doc.texts) if (inside(t)) this.sel.add(t.id);
  }

  /**
   * Bracket tool: the atoms in the box become one repeat unit with its own count. Polymer * ends
   * caught in the box stay outside the brackets, as chemists draw *–[unit]ₙ–*.
   */
  private bracketBox(a: Point, b: Point) {
    const x0 = Math.min(a.x, b.x), x1 = Math.max(a.x, b.x), y0 = Math.min(a.y, b.y), y1 = Math.max(a.y, b.y);
    const ids = new Set(this.doc.atoms.filter((at) => at.x >= x0 && at.x <= x1 && at.y >= y0 && at.y <= y1 &&
      !(at.el === '*' && bondsOf(this.doc, at.id).length <= 1)).map((at) => at.id));
    this.addBracket(ids);
  }

  /**
   * Bracket tool, click mode: the first bond clicked is remembered; the second one closes the unit,
   * which is everything between the two bonds (with its side groups).
   */
  private bracketBondClick(bondId: number) {
    const first = this.pendingBracketBond;
    if (first === null || !this.doc.bonds.some((b) => b.id === first)) {
      this.pendingBracketBond = bondId;
      this.render();
      this.onHover({ title: 'Bracket', detail: 'Now click the bond on the other side of the repeat unit. Click the same bond again or press Esc to cancel.' });
      return;
    }
    this.pendingBracketBond = null;
    if (first === bondId) { this.render(); this.emitHover(); return; }
    const ids = this.unitBetween(first, bondId);
    if (!ids) {
      this.render();
      this.onWarn('These two bonds do not cut out one unit: the atoms between them are still joined to the rest another way (for example through a ring). Click two chain bonds, one on each side of the unit.');
      return;
    }
    this.addBracket(ids);
  }

  /** Atoms cut off by removing bonds b1 and b2 that touch both of them, or null if there is no such part. */
  private unitBetween(b1: number, b2: number): Set<number> | null {
    const B1 = this.doc.bonds.find((b) => b.id === b1)!, B2 = this.doc.bonds.find((b) => b.id === b2)!;
    const adj = new Map<number, number[]>(this.doc.atoms.map((a) => [a.id, []]));
    for (const b of this.doc.bonds) if (b.id !== b1 && b.id !== b2) { adj.get(b.a)!.push(b.b); adj.get(b.b)!.push(b.a); }
    const part = (start: number) => {
      const seen = new Set([start]);
      const stack = [start];
      while (stack.length) for (const n of adj.get(stack.pop()!)!) if (!seen.has(n)) { seen.add(n); stack.push(n); }
      return seen;
    };
    let best: Set<number> | null = null;
    for (const start of [B1.a, B1.b]) {
      const c = part(start);
      const ends1 = +c.has(B1.a) + +c.has(B1.b), ends2 = +c.has(B2.a) + +c.has(B2.b);
      if (ends1 === 1 && ends2 === 1 && (!best || c.size < best.size)) best = c;
    }
    return best;
  }

  private addBracket(ids: Set<number>) {
    const problem = bracketProblem(this.doc, ids);
    if (problem) { this.render(); this.onWarn(problem); return; }
    if ((this.doc.brackets ?? []).some((g) => g.atoms.some((id) => ids.has(id)))) {
      this.render();
      this.onWarn('Some of these atoms are already inside a bracketed unit. Each atom can belong to one repeat unit: click that bracket to change its count, or erase it first.');
      return;
    }
    // the next unused letter, so blocks of a copolymer get n, m, p… unless the user types numbers
    const used = new Set((this.doc.brackets ?? []).map((g) => g.sub));
    for (const s of this.doc.atoms) if (s.el === '*' && repeatUnitEnds(this.doc).has(s.id)) used.add(s.sub ?? 'n');
    const sub = ['n', 'm', 'p', 'q', 'x', 'y', 'z'].find((l) => !used.has(l)) ?? 'n';
    let id = 0;
    this.commit((d) => { id = d.nextId++; (d.brackets ??= []).push({ id, atoms: [...ids], sub }); });
    if (id) this.editBracket(id);
  }

  private pendingBracketBond: number | null = null;

  /** Type the count for one bracketed repeat unit (n, m, 50, x…). */
  editBracket(id: number) {
    const m = bracketMarks(this.doc).find((x) => x.group === id && x.sub);
    const g = (this.doc.brackets ?? []).find((x) => x.id === id);
    if (!g) return;
    const at = m?.sub ? { x: m.sub.x, y: m.sub.y - 2 } : centroid(this.doc.atoms.filter((a) => g.atoms.includes(a.id)));
    this.openInput(at, g.sub, (v) => this.commit((d) => { const x = (d.brackets ?? []).find((q) => q.id === id); if (x) x.sub = v.trim() || 'n'; }), 'n');
  }

  /** Set the count of every repeat unit at once (homopolymer n), or of one bracket group. */
  setBracketCount(id: number, n: string) {
    this.commit((d) => { const g = (d.brackets ?? []).find((q) => q.id === id); if (g) g.sub = n || 'n'; });
  }

  /** Dropping a single dragged atom onto another atom merges them. */
  private mergeDroppedAtom(d: Doc, drag: Drag) {
    if (!drag.moveIds || drag.moveIds.size !== 1) return;
    const [id] = [...drag.moveIds];
    const moved = atomById(d, id);
    if (!moved) return;
    const target = findAtomNear(d, moved, BOND * 0.3, new Set([id]));
    if (!target) return;
    for (const b of bondsOf(d, id)) {
      const other = b.a === id ? b.b : b.a;
      if (other === target.id) continue;
      if (!d.bonds.some((x) => (x.a === other && x.b === target.id) || (x.b === other && x.a === target.id))) {
        addBond(d, target.id, other, b.order, b.stereo);
      }
    }
    removeObjects(d, new Set([id]));
  }

  selectAll() {
    for (const o of [...this.doc.atoms, ...this.doc.bonds, ...this.doc.arrows, ...this.doc.texts]) this.sel.add(o.id);
    this.render();
    this.onChange();
  }

  deleteSelection() {
    if (!this.sel.size) return;
    const ids = new Set(this.sel);
    this.commit((d) => {
      const touched = d.bonds.filter((b) => ids.has(b.id)).flatMap((b) => [b.a, b.b]);
      removeObjects(d, ids);
      pruneLoneCarbons(d, touched);
    });
    this.sel.clear();
  }

  private eraseAt(p: Point, first: boolean) {
    const hit = this.hitTest(p);
    if (!hit) { if (first) this.render(); return; }
    const touched = hit.type === 'bond' ? (() => { const b = this.doc.bonds.find((x) => x.id === hit.id)!; return [b.a, b.b]; })() : [];
    removeObjects(this.doc, new Set([hit.id]));
    pruneLoneCarbons(this.doc, touched);
    this.version++;
    this.hover = null;
    this.render();
  }

  /** Atom ids the transform actions should act on: selection if any, otherwise everything. */
  private targetIds(): { atoms: Set<number>; objs: Set<number> } {
    const atoms = this.movableAtoms();
    const objs = new Set([...this.sel].filter((id) => this.doc.arrows.some((a) => a.id === id) || this.doc.texts.some((t) => t.id === id)));
    if (!atoms.size && !objs.size) {
      return { atoms: new Set(this.doc.atoms.map((a) => a.id)), objs: new Set([...this.doc.arrows, ...this.doc.texts].map((o) => o.id)) };
    }
    return { atoms, objs };
  }

  rotate(deg: number) {
    const { atoms, objs } = this.targetIds();
    const c = centroid(this.doc.atoms.filter((a) => atoms.has(a.id)));
    const ang = (deg * Math.PI) / 180;
    this.commit((d) => {
      for (const a of d.atoms) if (atoms.has(a.id)) Object.assign(a, rotatePoint(a, c, ang));
      for (const ar of d.arrows) if (objs.has(ar.id)) {
        const p1 = rotatePoint({ x: ar.x1, y: ar.y1 }, c, ang), p2 = rotatePoint({ x: ar.x2, y: ar.y2 }, c, ang);
        Object.assign(ar, { x1: p1.x, y1: p1.y, x2: p2.x, y2: p2.y });
      }
    });
  }

  flip(axis: 'h' | 'v') {
    const { atoms } = this.targetIds();
    const c = centroid(this.doc.atoms.filter((a) => atoms.has(a.id)));
    this.commit((d) => {
      for (const a of d.atoms) if (atoms.has(a.id)) {
        if (axis === 'h') a.x = 2 * c.x - a.x;
        else a.y = 2 * c.y - a.y;
      }
    });
  }

  /** Re-layout structures with OpenChemLib's coordinate generator. */
  cleanUp() {
    const { atoms } = this.targetIds();
    if (!atoms.size) return;
    const part = subDoc(this.doc, atoms);
    const coords = cleanCoordinates(part);
    this.commit((d) => { for (const a of d.atoms) { const c = coords.get(a.id); if (c) { a.x = c.x; a.y = c.y; } } });
  }

  /**
   * Turn vinyl monomers into repeat units: the C=C becomes a single bond and both carbons
   * get a polymer end ('*'), drawn with brackets.
   */
  polymerize(): boolean {
    const { atoms } = this.targetIds();
    const an = this.analysis();
    const cand = this.doc.bonds.filter((b) => b.order === 2 && atoms.has(b.a) && atoms.has(b.b) && !an.ringCenter.has(b.id) &&
      atomById(this.doc, b.a)!.el === 'C' && atomById(this.doc, b.b)!.el === 'C');
    if (!cand.length) return false;
    // Prefer terminal vinyl groups (CH2=)
    cand.sort((x, y) => this.vinylScore(y) - this.vinylScore(x));
    const target = cand[0];
    this.commit((d) => {
      const b = d.bonds.find((x) => x.id === target.id)!;
      b.order = 1;
      for (const id of [b.a, b.b]) {
        const v = dirVec(freeDirection(d, id));
        const A = atomById(d, id)!;
        const s = addAtom(d, '*', A.x + v.x * BOND, A.y + v.y * BOND);
        addBond(d, id, s.id);
      }
    });
    this.sel.clear();
    this.render();
    this.onChange();
    return true;
  }

  private vinylScore(b: { a: number; b: number }) {
    return 4 - neighbors(this.doc, b.a).length - neighbors(this.doc, b.b).length;
  }

  setRepeatCount(n: string) {
    this.commit((d) => { for (const a of d.atoms) if (a.el === '*') a.sub = n || 'n'; });
  }

  // ---------- text editing ----------

  private openInput(world: Point, value: string, onCommit: (v: string) => void, placeholder = '') {
    this.closeInput();
    const inp = document.createElement('input');
    inp.className = 'canvas-input';
    inp.value = value;
    inp.placeholder = placeholder;
    inp.style.left = `${world.x * this.zoom + this.tx}px`;
    inp.style.top = `${world.y * this.zoom + this.ty - 14}px`;
    inp.addEventListener('keydown', (e) => {
      e.stopPropagation();
      if (e.key === 'Enter') this.commitInput();
      if (e.key === 'Escape') this.closeInput();
    });
    (inp as HTMLInputElement & { _commit?: (v: string) => void })._commit = onCommit;
    this.host.appendChild(inp);
    this.input = inp;
    inp.focus();
    inp.select();
  }

  private commitInput() {
    const inp = this.input as (HTMLInputElement & { _commit?: (v: string) => void }) | null;
    if (!inp) return;
    const v = inp.value;
    this.closeInput();
    inp._commit?.(v);
  }

  private closeInput() {
    this.input?.remove();
    this.input = null;
  }

  private newText(p: Point) {
    this.openInput(p, '', (v) => {
      if (!v.trim()) return;
      this.commit((d) => d.texts.push({ id: d.nextId++, x: p.x, y: p.y + 5, text: v }));
    }, 'Text…');
  }

  editText(id: number) {
    const t = this.doc.texts.find((x) => x.id === id) as TextItem;
    this.openInput({ x: t.x, y: t.y - 5 }, t.text, (v) => {
      this.commit((d) => {
        if (!v.trim()) d.texts = d.texts.filter((x) => x.id !== id);
        else d.texts.find((x) => x.id === id)!.text = v;
      });
    });
  }

  /** Text tool on an atom: type an element symbol, or n for a polymer end. */
  editAtom(id: number) {
    const a = atomById(this.doc, id)!;
    if (a.el === '*') {
      this.openInput(a, a.sub ?? 'n', (v) => this.commit((d) => { atomById(d, id)!.sub = v.trim() || 'n'; }));
      return;
    }
    this.openInput(a, a.el, (v) => this.setAtomLabel(id, v.trim()), 'Element');
  }

  /** Set an atom's element; also accepts a group abbreviation (OMe, CO2H, Ph…) via `expandAbbreviation`. */
  setAtomLabel(id: number, label: string) {
    if (!label) return;
    const el = label.length <= 2 ? label[0].toUpperCase() + label.slice(1).toLowerCase() : label;
    if (this.isElement(el) || el === '*') {
      this.commit((d) => { const a = atomById(d, id)!; a.el = el; if (el === '*') delete a.radical; });
      return;
    }
    if (!this.expandAbbreviation(id, label)) this.onStatus(`Unknown element or group “${label}”.`);
  }

  isElement: (el: string) => boolean = () => true;
  abbreviation: (label: string) => Doc | null = () => null;

  private expandAbbreviation(id: number, label: string): boolean {
    const frag = this.abbreviation(label);
    if (!frag) return false;
    const nb = neighbors(this.doc, id);
    const A = atomById(this.doc, id)!;
    if (nb.length > 1) return false;
    this.commit((d) => {
      if (nb.length === 0) {
        removeObjects(d, new Set([id]));
        const pl = placeTemplate(d, frag, 'group', A);
        mergeInto(d, pl.frag, 0, 0, pl.map);
      } else {
        const N = nb[0];
        const ang = Math.atan2(A.y - N.y, A.x - N.x);
        removeObjects(d, new Set([id]));
        const pl = placeTemplate(d, frag, 'group', A, N.id, ang);
        mergeInto(d, pl.frag, 0, 0, pl.map);
      }
    });
    return true;
  }

  // ---------- keyboard ----------

  private keyDown(e: KeyboardEvent) {
    const target = e.target as HTMLElement;
    if (target && (target.tagName === 'INPUT' || target.tagName === 'TEXTAREA' || target.isContentEditable)) return;
    const mod = e.metaKey || e.ctrlKey;
    if (e.code === 'Space' && !mod) { this.spaceDown = true; this.svg.style.cursor = 'grab'; e.preventDefault(); return; }
    if (mod && e.key.toLowerCase() === 'z') { e.preventDefault(); if (e.shiftKey) this.redo(); else this.undo(); return; }
    if (mod && e.key.toLowerCase() === 'y') { e.preventDefault(); this.redo(); return; }
    if (mod && e.key.toLowerCase() === 'a') { e.preventDefault(); this.selectAll(); return; }
    if (mod && (e.key === '=' || e.key === '+')) { e.preventDefault(); this.zoomBy(1.2); return; }
    if (mod && e.key === '-') { e.preventDefault(); this.zoomBy(1 / 1.2); return; }
    if (mod && e.key === '0') { e.preventDefault(); this.fitView(); return; }
    if (mod) return;
    if (e.key === 'Escape') {
      this.sel.clear();
      this.closeInput();
      this.setTool({ kind: 'select' });
      this.onChange();
      return;
    }
    if (e.key === 'Delete' || e.key === 'Backspace') {
      e.preventDefault();
      if (this.sel.size) this.deleteSelection();
      else if (this.hover) {
        const h = this.hover;
        this.commit((d) => {
          const touched = h.type === 'bond' ? (() => { const b = d.bonds.find((x) => x.id === h.id)!; return [b.a, b.b]; })() : [];
          removeObjects(d, new Set([h.id]));
          pruneLoneCarbons(d, touched);
        });
        this.hover = null;
      }
      return;
    }
    const h = this.hover;
    if (h?.type === 'atom') {
      const el = HOTKEY_ELEMENTS[e.key];
      if (el) { e.preventDefault(); this.commit((d) => { const a = atomById(d, h.id)!; a.el = el; if (el === '*') delete a.radical; }); return; }
      if (e.key === '.') { e.preventDefault(); this.cycleRadical(h.id); return; }
      if (e.key === '+' || e.key === '=') { this.commit((d) => { atomById(d, h.id)!.charge += 1; }); return; }
      if (e.key === '-') { this.commit((d) => { atomById(d, h.id)!.charge -= 1; }); return; }
      if (e.key === 'Enter') { this.editAtom(h.id); return; }
    }
    if (h?.type === 'bond' && ['1', '2', '3'].includes(e.key)) {
      const order = parseInt(e.key, 10) as 1 | 2 | 3;
      this.commit((d) => { const b = d.bonds.find((x) => x.id === h.id)!; b.order = order; b.stereo = 'none'; });
      return;
    }
    if (h && (h.type === 'atom' ? HOTKEY_ELEMENTS[e.key] : ['1', '2', '3'].includes(e.key))) return;
    const shortcuts: Record<string, Tool> = {
      v: { kind: 'select' }, e: { kind: 'erase' }, b: { kind: 'bond', order: 1, stereo: 'none' }, g: { kind: 'chain' },
      r: { kind: 'ring', size: 6, aromatic: true }, t: { kind: 'text' }, a: { kind: 'arrow', arrow: 'forward' }, '[': { kind: 'bracket' },
      '1': { kind: 'bond', order: 1, stereo: 'none' }, '2': { kind: 'bond', order: 2, stereo: 'none' }, '3': { kind: 'bond', order: 3, stereo: 'none' },
    };
    const tool = shortcuts[e.key];
    if (tool) this.setTool(tool);
  }

  // ---------- clipboard ----------

  selectionDoc(): Doc {
    const atoms = this.movableAtoms();
    return subDoc(this.doc, atoms, this.sel);
  }

  private copy(e: ClipboardEvent, cut: boolean) {
    const target = e.target as HTMLElement;
    if (target && (target.tagName === 'INPUT' || target.tagName === 'TEXTAREA')) return;
    if (!this.sel.size) return;
    const part = this.selectionDoc();
    const smiles = part.atoms.length ? toSmiles(part) : '';
    this.clipboard = { json: JSON.stringify(part), smiles };
    e.clipboardData?.setData('text/plain', smiles);
    e.preventDefault();
    if (cut) this.deleteSelection();
    this.onStatus(cut ? 'Cut to clipboard.' : `Copied${smiles ? ` (SMILES: ${smiles})` : ''}.`);
  }

  private paste(e: ClipboardEvent) {
    const target = e.target as HTMLElement;
    if (target && (target.tagName === 'INPUT' || target.tagName === 'TEXTAREA')) return;
    const text = e.clipboardData?.getData('text/plain') ?? '';
    e.preventDefault();
    if (this.clipboard && (text === this.clipboard.smiles || !text)) {
      const frag = JSON.parse(this.clipboard.json) as Doc;
      this.insertFragment(frag, true);
      return;
    }
    if (text.trim()) {
      try { this.insertText(text); }
      catch { this.onStatus('Could not read the clipboard as SMILES or a MOL file.'); }
    }
  }

  /** Insert SMILES, MOL or ChemDraw CDXML text at the view centre. Throws if it cannot be parsed. */
  insertText(text: string) {
    if (/<CDXML[\s>]/.test(text)) {
      const { doc, notes } = fromCdxml(text);
      this.insertFragment(doc, false);
      if (notes.length) this.onStatus(`ChemDraw import: ${notes.join(' ')}`);
      return;
    }
    const frag = /V2000|V3000|M {2}END/.test(text) ? fromMolfile(text) : fromSmiles(text.trim().split(/\s+/)[0]);
    orientPolymer(frag);
    this.insertFragment(frag, false);
  }

  insertFragment(frag: Doc, offset: boolean) {
    if (!frag.atoms.length && !frag.arrows.length && !frag.texts.length) return;
    const c = centroid(frag.atoms.length ? frag.atoms : frag.texts);
    const target = offset ? null : this.freeSpotFor(frag);
    const dx = offset ? BOND * 0.8 : target!.x - c.x, dy = offset ? BOND * 0.8 : target!.y - c.y;
    if (offset) {
      // move the clipboard copy so repeated pastes cascade
      for (const a of frag.atoms) { a.x += dx; a.y += dy; }
      for (const t of frag.texts) { t.x += dx; t.y += dy; }
      for (const a of frag.arrows) { a.x1 += dx; a.y1 += dy; a.x2 += dx; a.y2 += dy; }
      if (this.clipboard) this.clipboard.json = JSON.stringify(frag);
    }
    let map = new Map<number, number>();
    this.commit((d) => { map = mergeInto(d, frag, offset ? 0 : dx, offset ? 0 : dy); });
    this.sel = new Set(map.values());
    if (target) this.ensureVisible(target);
    this.render();
    this.onChange();
  }

  /** Where to drop new content: the view centre on an empty canvas, otherwise right of what is there. */
  private freeSpotFor(frag: Doc): Point {
    const vc = this.viewCenter();
    const b = bounds(this.doc);
    if (!b) return vc;
    const fb = bounds(frag)!;
    return { x: b.x + b.w + BOND * 3.5 + fb.w / 2, y: b.y + b.h / 2 };
  }

  private ensureVisible(p: Point) {
    const r = this.svg.getBoundingClientRect();
    const sx = p.x * this.zoom + this.tx, sy = p.y * this.zoom + this.ty;
    const m = 120;
    if (sx < m || sx > r.width - m || sy < m || sy > r.height - m) {
      this.tx = r.width / 2 - p.x * this.zoom;
      this.ty = r.height / 2 - p.y * this.zoom;
    }
  }
}
