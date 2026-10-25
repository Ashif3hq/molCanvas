// SVG renderer for MolCanvas documents. Produces markup strings so the same code
// drives the live canvas, library thumbnails and SVG/PNG export.
import { Analysis, analyze } from './chem';
import { Atom, Bond, BOND, Doc, Point, Arrow, TextItem, fragments, repeatUnitEnds } from './model';

export interface Theme {
  fg: string;
  bg: string;
  select: string;
  hover: string;
  error: string;
}

export const LIGHT: Theme = { fg: '#1d1d1f', bg: '#ffffff', select: '#0a84ff', hover: '#0a84ff', error: '#ff3b30' };
export const DARK: Theme = { fg: '#f2f2f7', bg: '#1c1c1e', select: '#0a84ff', hover: '#64d2ff', error: '#ff453a' };

export interface RenderOpts {
  theme: Theme;
  colorAtoms: boolean;
  selection?: Set<number>;
  hover?: number | null;
  analysis?: Analysis;
  showErrors?: boolean;
}

export const FONT = 16;
const LINE = 1.8;
const DOUBLE_GAP = 6.5;
const FONT_FAMILY = "'Helvetica Neue', Helvetica, Arial, sans-serif";

const ELEMENT_COLORS: Record<string, string> = {
  N: '#2f5bea', O: '#e0322b', S: '#c49a00', P: '#e07000', F: '#1e9e3a', Cl: '#1e9e3a', Br: '#a52a2a', I: '#8a1fa0',
  B: '#c9733a', Si: '#9a7752', Se: '#c47a00', Li: '#8a4fd6', Na: '#8a4fd6', K: '#8a4fd6', Mg: '#3c9b3c',
};

const esc = (s: string) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const f = (n: number) => (Math.round(n * 100) / 100).toString();

export function atomColor(el: string, opts: RenderOpts): string {
  if (!opts.colorAtoms) return opts.theme.fg;
  return ELEMENT_COLORS[el] ?? opts.theme.fg;
}

function charWidth(ch: string): number {
  if (/[0-9]/.test(ch)) return 0.56;
  if (ch === 'I' || ch === 'l' || ch === 'i') return 0.3;
  if (/[A-Z]/.test(ch)) return 0.7;
  return 0.55;
}

export function textWidth(s: string, size = FONT): number {
  let w = 0;
  for (const ch of s) w += charWidth(ch);
  return w * size;
}

interface LabelInfo {
  shown: boolean;
  symbol: string;
  hCount: number;
  hSide: 'left' | 'right';
  halfW: number; // half width of the element symbol
}

function chargeText(c: number): string {
  if (!c) return '';
  const mag = Math.abs(c);
  return `${mag > 1 ? mag : ''}${c > 0 ? '+' : '−'}`;
}

export function labelInfo(doc: Doc, atom: Atom, nbrs: Atom[], an: Analysis, bracketed = false): LabelInfo {
  const h = an.hCount.get(atom.id) ?? 0;
  let shown = atom.el !== 'C' || nbrs.length === 0 || !!atom.iso;
  if (atom.el === '*') shown = !bracketed;
  let sx = 0;
  for (const n of nbrs) sx += n.x - atom.x;
  const hSide: 'left' | 'right' = nbrs.length && sx > 1 ? 'left' : 'right';
  const symbol = (atom.iso ? '' : '') + atom.el;
  // 'R' marks a generic substituent (used for attachment points in group previews)
  return { shown, symbol, hCount: shown && atom.el !== '*' && atom.el !== 'R' ? h : 0, hSide, halfW: textWidth(symbol) / 2 };
}

/** How far to pull a bond end back from a labelled atom, given the bond direction (ux, uy). */
function clipRadius(li: LabelInfo, ux: number, uy: number): number {
  if (!li.shown) return 0;
  let hw = li.halfW + 2.5;
  const hh = FONT * 0.55;
  // Hydrogens printed on the side the bond leaves from need extra room.
  if (li.hCount && ((li.hSide === 'right' && ux > 0.3) || (li.hSide === 'left' && ux < -0.3))) hw += textWidth('H2');
  const r = 1 / Math.sqrt((ux / hw) ** 2 + (uy / hh) ** 2);
  return Math.min(r, BOND * 0.4);
}

interface Ctx {
  doc: Doc;
  opts: RenderOpts;
  an: Analysis;
  byId: Map<number, Atom>;
  nbrs: Map<number, Atom[]>;
  labels: Map<number, LabelInfo>;
}

function buildCtx(doc: Doc, opts: RenderOpts): Ctx {
  const an = opts.analysis ?? analyze(doc);
  const byId = new Map<number, Atom>();
  const nbrs = new Map<number, Atom[]>();
  for (const a of doc.atoms) { byId.set(a.id, a); nbrs.set(a.id, []); }
  for (const b of doc.bonds) {
    const A = byId.get(b.a), B = byId.get(b.b);
    if (!A || !B) continue;
    nbrs.get(A.id)!.push(B);
    nbrs.get(B.id)!.push(A);
  }
  const labels = new Map<number, LabelInfo>();
  const ends = repeatUnitEnds(doc);
  for (const a of doc.atoms) labels.set(a.id, labelInfo(doc, a, nbrs.get(a.id)!, an, ends.has(a.id)));
  return { doc, opts, an, byId, nbrs, labels };
}

function line(x1: number, y1: number, x2: number, y2: number, color: string, w = LINE, extra = ''): string {
  return `<line x1="${f(x1)}" y1="${f(y1)}" x2="${f(x2)}" y2="${f(y2)}" stroke="${color}" stroke-width="${w}" stroke-linecap="round"${extra}/>`;
}

function bondSvg(ctx: Ctx, b: Bond): string {
  const A = ctx.byId.get(b.a), B = ctx.byId.get(b.b);
  if (!A || !B) return '';
  const dx = B.x - A.x, dy = B.y - A.y;
  const len = Math.hypot(dx, dy) || 1;
  const ux = dx / len, uy = dy / len;
  const px = -uy, py = ux;
  const ra = clipRadius(ctx.labels.get(A.id)!, ux, uy);
  const rb = clipRadius(ctx.labels.get(B.id)!, -ux, -uy);
  const x1 = A.x + ux * ra, y1 = A.y + uy * ra;
  const x2 = B.x - ux * rb, y2 = B.y - uy * rb;
  const fg = ctx.opts.theme.fg;
  let out = '';

  if (b.order === 1) {
    if (b.stereo === 'wedge') {
      const w = 4.2;
      out += `<polygon points="${f(x1)},${f(y1)} ${f(x2 + px * w)},${f(y2 + py * w)} ${f(x2 - px * w)},${f(y2 - py * w)}" fill="${fg}" stroke="${fg}" stroke-width="0.6" stroke-linejoin="round"/>`;
    } else if (b.stereo === 'hash') {
      const n = Math.max(4, Math.round(len / 5.5));
      for (let i = 0; i <= n; i++) {
        const t = i / n;
        const cx = x1 + (x2 - x1) * t, cy = y1 + (y2 - y1) * t;
        const w = 0.6 + 3.8 * t;
        out += line(cx + px * w, cy + py * w, cx - px * w, cy - py * w, fg, 1.4);
      }
    } else if (b.stereo === 'wavy') {
      const seg = Math.max(4, Math.round(len / 6));
      let d = `M${f(x1)},${f(y1)}`;
      for (let i = 1; i <= seg; i++) {
        const t0 = (i - 0.5) / seg, t1 = i / seg;
        const s = i % 2 ? 1 : -1;
        const cx = x1 + (x2 - x1) * t0 + px * 4 * s, cy = y1 + (y2 - y1) * t0 + py * 4 * s;
        d += ` Q${f(cx)},${f(cy)} ${f(x1 + (x2 - x1) * t1)},${f(y1 + (y2 - y1) * t1)}`;
      }
      out += `<path d="${d}" fill="none" stroke="${fg}" stroke-width="${LINE}" stroke-linecap="round"/>`;
    } else {
      out += line(x1, y1, x2, y2, fg);
    }
    return out;
  }

  if (b.order === 3) {
    const g = DOUBLE_GAP * 0.95;
    out += line(x1, y1, x2, y2, fg);
    out += line(x1 + px * g, y1 + py * g, x2 + px * g, y2 + py * g, fg);
    out += line(x1 - px * g, y1 - py * g, x2 - px * g, y2 - py * g, fg);
    return out;
  }

  // Double bond of unknown E/Z geometry: drawn as two crossed lines.
  if (b.stereo === 'wavy') {
    const g = DOUBLE_GAP / 2;
    out += line(x1 + px * g, y1 + py * g, x2 - px * g, y2 - py * g, fg);
    out += line(x1 - px * g, y1 - py * g, x2 + px * g, y2 + py * g, fg);
    return out;
  }

  // Double bond: decide whether to centre it or offset the second line to one side.
  const side = doubleBondSide(ctx, b, A, B, px, py);
  if (side === 0) {
    const g = DOUBLE_GAP / 2;
    out += line(x1 + px * g, y1 + py * g, x2 + px * g, y2 + py * g, fg);
    out += line(x1 - px * g, y1 - py * g, x2 - px * g, y2 - py * g, fg);
  } else {
    out += line(x1, y1, x2, y2, fg);
    const g = DOUBLE_GAP * side;
    const la = ctx.labels.get(A.id)!.shown, lb = ctx.labels.get(B.id)!.shown;
    const sa = la ? 0 : len * 0.13, sb = lb ? 0 : len * 0.13;
    out += line(x1 + ux * sa + px * g, y1 + uy * sa + py * g, x2 - ux * sb + px * g, y2 - uy * sb + py * g, fg);
  }
  return out;
}

/** -1, 0 (centred) or +1: which side of the bond (relative to the perpendicular) the second line goes. */
function doubleBondSide(ctx: Ctx, b: Bond, A: Atom, B: Atom, px: number, py: number): number {
  const ring = ctx.an.ringCenter.get(b.id);
  if (ring) {
    const s = (ring.x - A.x) * px + (ring.y - A.y) * py;
    return s >= 0 ? 1 : -1;
  }
  const na = ctx.nbrs.get(A.id)!.filter((n) => n.id !== B.id);
  const nb = ctx.nbrs.get(B.id)!.filter((n) => n.id !== A.id);
  const la = ctx.labels.get(A.id)!.shown, lb = ctx.labels.get(B.id)!.shown;
  if ((na.length === 0 && (nb.length !== 1 || lb)) || (nb.length === 0 && (na.length !== 1 || la))) return 0;
  if (na.length === 0 && nb.length === 0) return 0;
  let score = 0;
  for (const n of na) score += Math.sign((n.x - A.x) * px + (n.y - A.y) * py);
  for (const n of nb) score += Math.sign((n.x - B.x) * px + (n.y - B.y) * py);
  if (score === 0) return 0;
  return score > 0 ? 1 : -1;
}

function atomSvg(ctx: Ctx, a: Atom): string {
  const li = ctx.labels.get(a.id)!;
  const color = atomColor(a.el, ctx.opts);
  let out = '';
  const baseY = a.y + FONT * 0.36;
  if (li.shown) {
    const label = a.el;
    if (a.iso) {
      out += `<text x="${f(a.x - li.halfW - 1)}" y="${f(a.y - FONT * 0.2)}" font-size="${FONT * 0.65}" text-anchor="end" fill="${color}">${a.iso}</text>`;
    }
    out += `<text x="${f(a.x)}" y="${f(baseY)}" text-anchor="middle" fill="${color}">${esc(label)}</text>`;
    if (li.hCount > 0) {
      const sub = li.hCount > 1 ? `<tspan font-size="${FONT * 0.72}" dy="${FONT * 0.25}">${li.hCount}</tspan>` : '';
      if (li.hSide === 'right') {
        out += `<text x="${f(a.x + li.halfW)}" y="${f(baseY)}" text-anchor="start" fill="${color}">H${sub}</text>`;
      } else {
        out += `<text x="${f(a.x - li.halfW)}" y="${f(baseY)}" text-anchor="end" fill="${color}">H${sub}</text>`;
      }
    }
  }
  if (a.charge) {
    let cx = a.x + (li.shown ? li.halfW + 1 : 5);
    if (li.shown && li.hCount && li.hSide === 'right') cx += textWidth(li.hCount > 1 ? 'H2' : 'H');
    const cy = a.y - (li.shown ? FONT * 0.38 : 7);
    out += `<text x="${f(cx)}" y="${f(cy)}" font-size="${FONT * 0.75}" text-anchor="start" fill="${color}">${chargeText(a.charge)}</text>`;
  }
  if (ctx.opts.showErrors && ctx.an.badValence.has(a.id)) {
    out += `<circle cx="${f(a.x)}" cy="${f(a.y)}" r="11" fill="none" stroke="${ctx.opts.theme.error}" stroke-width="1.4" stroke-dasharray="3 2"/>`;
  }
  return out;
}

/** Polymer brackets: drawn across every bond to a '*' end atom, with the subscript on the right-most one. */
function bracketsSvg(ctx: Ctx): string {
  let out = '';
  const fg = ctx.opts.theme.fg;
  const frags = fragments(ctx.doc);
  const ends = repeatUnitEnds(ctx.doc);
  for (const comp of frags) {
    const stars = comp.map((id) => ctx.byId.get(id)!).filter((a) => ends.has(a.id));
    if (!stars.length) continue;
    const last = stars.reduce((m, s) => (s.x > m.x ? s : m), stars[0]);
    for (const s of stars) {
      const X = ctx.nbrs.get(s.id)![0];
      const dx = X.x - s.x, dy = X.y - s.y;
      const len = Math.hypot(dx, dy) || 1;
      const ux = dx / len, uy = dy / len;
      const px = -uy, py = ux;
      const t = 0.42;
      const cx = s.x + dx * t, cy = s.y + dy * t;
      const h = BOND * 0.5, serif = 5;
      const e1 = { x: cx + px * h, y: cy + py * h }, e2 = { x: cx - px * h, y: cy - py * h };
      out += `<path d="M${f(e1.x + ux * serif)},${f(e1.y + uy * serif)} L${f(e1.x)},${f(e1.y)} L${f(e2.x)},${f(e2.y)} L${f(e2.x + ux * serif)},${f(e2.y + uy * serif)}" fill="none" stroke="${fg}" stroke-width="${LINE}" stroke-linejoin="miter"/>`;
      if (s === last) {
        const low = e1.y > e2.y ? e1 : e2;
        const sub = s.sub ?? 'n';
        out += `<text x="${f(low.x + 4)}" y="${f(low.y + 11)}" font-size="${FONT * 0.85}" font-style="italic" fill="${fg}" data-sub="${s.id}">${esc(sub)}</text>`;
      }
    }
  }
  return out;
}

function arrowSvg(a: Arrow, color: string): string {
  const dx = a.x2 - a.x1, dy = a.y2 - a.y1;
  const len = Math.hypot(dx, dy) || 1;
  const ux = dx / len, uy = dy / len;
  const px = -uy, py = ux;
  const head = 11, hw = 4.5;
  const tri = (tx: number, ty: number, dirx: number, diry: number) => {
    const bx = tx - dirx * head, by = ty - diry * head;
    return `<polygon points="${f(tx)},${f(ty)} ${f(bx - diry * hw)},${f(by + dirx * hw)} ${f(bx + diry * hw)},${f(by - dirx * hw)}" fill="${color}"/>`;
  };
  let out = '';
  if (a.kind === 'forward' || a.kind === 'resonance') {
    const s = a.kind === 'resonance' ? head * 0.8 : 0;
    out += line(a.x1 + ux * s, a.y1 + uy * s, a.x2 - ux * head * 0.8, a.y2 - uy * head * 0.8, color);
    out += tri(a.x2, a.y2, ux, uy);
    if (a.kind === 'resonance') out += tri(a.x1, a.y1, -ux, -uy);
  } else if (a.kind === 'equilibrium') {
    const g = 3.2;
    out += line(a.x1 - px * g, a.y1 - py * g, a.x2 - px * g, a.y2 - py * g, color);
    out += line(a.x2 - px * g, a.y2 - py * g, a.x2 - ux * head - px * (g + hw + 1), a.y2 - uy * head - py * (g + hw + 1), color);
    out += line(a.x1 + px * g, a.y1 + py * g, a.x2 + px * g, a.y2 + py * g, color);
    out += line(a.x1 + px * g, a.y1 + py * g, a.x1 + ux * head + px * (g + hw + 1), a.y1 + uy * head + py * (g + hw + 1), color);
  } else {
    const g = 3.5;
    const ex = a.x2 - ux * 7, ey = a.y2 - uy * 7;
    out += line(a.x1 + px * g, a.y1 + py * g, ex + px * g, ey + py * g, color);
    out += line(a.x1 - px * g, a.y1 - py * g, ex - px * g, ey - py * g, color);
    out += `<path d="M${f(a.x2 - ux * 13 + px * 9)},${f(a.y2 - uy * 13 + py * 9)} L${f(a.x2)},${f(a.y2)} L${f(a.x2 - ux * 13 - px * 9)},${f(a.y2 - uy * 13 - py * 9)}" fill="none" stroke="${color}" stroke-width="${LINE}" stroke-linejoin="miter"/>`;
  }
  return out;
}

/** Text with chemistry-aware subscripts: digits directly after a letter or ')' are subscripted (CH2Cl2). */
export function chemTextSpans(text: string): string {
  let out = '';
  const re = /([A-Za-z)\]])(\d+)/g;
  let last = 0;
  let m: RegExpExecArray | null;
  while ((m = re.exec(text))) {
    out += esc(text.slice(last, m.index + 1));
    out += `<tspan font-size="${FONT * 0.68}" dy="${FONT * 0.22}">${m[2]}</tspan><tspan dy="${-FONT * 0.22}">​</tspan>`;
    last = m.index + m[0].length;
  }
  out += esc(text.slice(last));
  return out;
}

function textSvg(t: TextItem, color: string): string {
  return `<text x="${f(t.x)}" y="${f(t.y)}" fill="${color}" font-size="${FONT * 0.9}">${chemTextSpans(t.text)}</text>`;
}

export function textBox(t: TextItem) {
  const w = textWidth(t.text, FONT * 0.9) * 0.9 + 4;
  return { x: t.x - 2, y: t.y - FONT * 0.9, w, h: FONT * 1.25 };
}

function highlights(ctx: Ctx): string {
  const sel = ctx.opts.selection ?? new Set<number>();
  const hover = ctx.opts.hover;
  const { theme } = ctx.opts;
  let out = '';
  for (const b of ctx.doc.bonds) {
    const A = ctx.byId.get(b.a)!, B = ctx.byId.get(b.b)!;
    const selected = sel.has(b.id) || (sel.has(b.a) && sel.has(b.b));
    if (selected) out += line(A.x, A.y, B.x, B.y, theme.select, 10, ' stroke-opacity="0.28"');
    if (hover === b.id) out += line(A.x, A.y, B.x, B.y, theme.hover, 10, ' stroke-opacity="0.35"');
  }
  for (const a of ctx.doc.atoms) {
    if (sel.has(a.id)) out += `<circle cx="${f(a.x)}" cy="${f(a.y)}" r="9" fill="${theme.select}" fill-opacity="0.3"/>`;
    if (hover === a.id) out += `<circle cx="${f(a.x)}" cy="${f(a.y)}" r="10" fill="${theme.hover}" fill-opacity="0.35"/>`;
  }
  for (const ar of ctx.doc.arrows) {
    if (sel.has(ar.id) || hover === ar.id) out += line(ar.x1, ar.y1, ar.x2, ar.y2, sel.has(ar.id) ? theme.select : theme.hover, 12, ' stroke-opacity="0.25"');
  }
  for (const t of ctx.doc.texts) {
    if (sel.has(t.id) || hover === t.id) {
      const bx = textBox(t);
      out += `<rect x="${f(bx.x)}" y="${f(bx.y)}" width="${f(bx.w)}" height="${f(bx.h)}" rx="3" fill="${sel.has(t.id) ? theme.select : theme.hover}" fill-opacity="0.2"/>`;
    }
  }
  return out;
}

/** Inner SVG markup for a document (no <svg> wrapper). */
export function renderDoc(doc: Doc, opts: RenderOpts): string {
  const ctx = buildCtx(doc, opts);
  let out = `<g class="hl">${highlights(ctx)}</g>`;
  out += `<g font-family="${FONT_FAMILY}" font-size="${FONT}">`;
  for (const b of doc.bonds) out += bondSvg(ctx, b);
  out += bracketsSvg(ctx);
  for (const a of doc.atoms) out += atomSvg(ctx, a);
  for (const ar of doc.arrows) out += arrowSvg(ar, opts.theme.fg);
  for (const t of doc.texts) out += textSvg(t, opts.theme.fg);
  out += '</g>';
  return out;
}

export function docBounds(doc: Doc, pad: number) {
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  const take = (p: Point, r = 0) => { x0 = Math.min(x0, p.x - r); y0 = Math.min(y0, p.y - r); x1 = Math.max(x1, p.x + r); y1 = Math.max(y1, p.y + r); };
  for (const a of doc.atoms) take(a, a.el === 'C' ? 4 : 12);
  for (const a of doc.arrows) { take({ x: a.x1, y: a.y1 }, 6); take({ x: a.x2, y: a.y2 }, 6); }
  for (const t of doc.texts) { const b = textBox(t); take({ x: b.x, y: b.y }); take({ x: b.x + b.w, y: b.y + b.h }); }
  if (x0 === Infinity) return { x: 0, y: 0, w: 100, h: 100 };
  return { x: x0 - pad, y: y0 - pad, w: x1 - x0 + 2 * pad, h: y1 - y0 + 2 * pad };
}

/** A complete standalone SVG document, used for export and thumbnails. */
export function standaloneSvg(doc: Doc, opts: RenderOpts & { background?: string | null; pad?: number; scale?: number }): string {
  const b = docBounds(doc, opts.pad ?? 16);
  const s = opts.scale ?? 1;
  const bg = opts.background ? `<rect x="${f(b.x)}" y="${f(b.y)}" width="${f(b.w)}" height="${f(b.h)}" fill="${opts.background}"/>` : '';
  return `<svg xmlns="http://www.w3.org/2000/svg" width="${f(b.w * s)}" height="${f(b.h * s)}" viewBox="${f(b.x)} ${f(b.y)} ${f(b.w)} ${f(b.h)}">${bg}${renderDoc(doc, { ...opts, selection: undefined, hover: null })}</svg>`;
}
