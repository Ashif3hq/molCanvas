// NMR panel in the inspector (predicted ¹H / ¹³C spectrum of the drawing, with a peak table) and the larger
// spectrum dialog with zoom, export and a searchable list of reference peaks.
import { Doc } from './model';
import { toOCL, fromSmiles } from './chem';
import { DARK, LIGHT, standaloneSvg } from './render';
import { saveFile } from './platform';
import {
  couplingText, formatShift, NMR_SOLVENTS, NmrResult, NmrSolvent, NmrSpectrum, Nucleus, predictNmr, sampleSpectrum, spectrumLines,
} from './nmr';
import { fmtRef, NMR_REFS, NmrRef, REF_CATEGORIES, RefCategory, RefValue, SOLVENT_PEAKS } from './nmr-ref';

export const PEAK_COLORS = ['#0a84ff', '#ff9f0a', '#30b158', '#bf5af2', '#ff375f', '#40c8e0', '#a2845e', '#5e5ce6'];
const STORE = 'molcanvas.nmr.v1';
const NUC_LABEL: Record<Nucleus, string> = { '1H': '¹H', '13C': '¹³C' };
const esc = (t: string) => t.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');

interface State {
  nucleus: Nucleus;
  solvent: NmrSolvent;
  mhz: number;
  open: boolean;
  pinned: string[];
  showSolvent: boolean;
  showLit: boolean;
  cat: RefCategory | 'all';
  query: string;
}
let st: State = { nucleus: '1H', solvent: 'CDCl3', mhz: 400, open: false, pinned: [], showSolvent: true, showLit: true, cat: 'all', query: '' };
try { const raw = localStorage.getItem(STORE); if (raw) st = { ...st, ...JSON.parse(raw) }; } catch { /* storage unavailable */ }
const persist = () => { try { localStorage.setItem(STORE, JSON.stringify(st)); } catch { /* ignore */ } };

interface Host {
  setMarks: (marks: Map<number, string> | null) => void;
  isDark: () => boolean;
  notify: (msg: string) => void;
}
let host: Host;
let doc: Doc | null = null;
let label = '';
const results: Partial<Record<Nucleus, NmrResult>> = {};
let match: NmrRef | null = null;
let dialog: HTMLDialogElement | null = null;
const zoom: Partial<Record<Nucleus, [number, number]>> = {};

export function initNmr(h: Host) {
  host = h;
  const btn = document.getElementById('nmr-toggle')!;
  const panel = document.getElementById('nmr')!;
  const sync = () => { btn.setAttribute('aria-expanded', String(st.open)); panel.hidden = !st.open; };
  sync();
  btn.addEventListener('click', () => { st.open = !st.open; persist(); sync(); renderPanel(); });
}

/** Called by the inspector with the drawing (or the selection) whenever it changes. */
export function updateNmr(d: Doc, name: string) {
  doc = d.atoms.some((a) => a.el !== '*') ? d : null;
  label = name;
  for (const k of Object.keys(results)) delete results[k as Nucleus];
  match = null;
  if (doc) match = literatureMatch(doc);
  if (st.open) renderPanel();
  if (dialog?.open) renderDialog();
}

function result(n: Nucleus): NmrResult {
  if (!doc) return { ok: false, reason: 'Nothing drawn.' };
  return (results[n] ??= predictNmr(doc, n, st.solvent));
}

// ---------- literature match ----------

let refCodes: Map<string, NmrRef> | null = null;
function literatureMatch(d: Doc): NmrRef | null {
  try {
    if (!refCodes) {
      refCodes = new Map();
      for (const r of NMR_REFS) if (r.smiles) refCodes.set(toOCL(fromSmiles(r.smiles)).mol.getIDCode(), r);
    }
    return refCodes.get(toOCL(d).mol.getIDCode()) ?? null;
  } catch { return null; }
}

// ---------- plotting ----------

interface Marker { v?: number; lo?: number; hi?: number; text: string; kind: 'solvent' | 'ref' | 'lit' }
interface Colors { fg: string; muted: string; grid: string; ref: string; lit: string; band: string; bg?: string }
const CSS_COLORS: Colors = { fg: 'var(--text)', muted: 'var(--muted)', grid: 'var(--border)', ref: '#bf5af2', lit: '#30b158', band: '#bf5af2' };

export function autoRange(sp: NmrSpectrum): [number, number] {
  const vals = sp.peaks.flatMap((p) => (p.range ? [p.range[0], p.range[1]] : [p.shift - p.width / 2, p.shift + p.width / 2]));
  if (sp.nucleus === '1H') {
    const hi = Math.max(5, Math.ceil(Math.max(...vals, 0) + 0.6));
    const lo = Math.min(0, Math.floor(Math.min(...vals, 0) - 0.3));
    return [hi, lo];
  }
  const hi = Math.max(80, Math.ceil((Math.max(...vals, 0) + 12) / 20) * 20);
  return [hi, Math.min(0, Math.floor((Math.min(...vals, 0) - 5) / 10) * 10)];
}

function tickStep(span: number): number {
  for (const s of [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50]) if (span / s <= 10) return s;
  return 100;
}

let plotSeq = 0;

export function plotSvg(sp: NmrSpectrum, o: { w: number; h: number; x0: number; x1: number; mhz: number; hl?: number | null; markers?: Marker[]; colors?: Colors; font?: number; markerText?: boolean }): string {
  const c = o.colors ?? CSS_COLORS;
  const clip = `nmr-clip-${++plotSeq}`;
  const markerText = o.markerText ?? true;
  const placedText: { a: number; b: number; y: number }[] = [];
  /** y for a marker label spanning [a, b] px, moved up a row while it would overlap another label. */
  const textY = (a: number, b: number, base: number) => {
    let y = base;
    while (placedText.some((p) => a < p.b + 4 && b > p.a - 4 && p.y === y)) y -= 11;
    placedText.push({ a, b, y });
    return y;
  };
  const fs = o.font ?? 10;
  const top = 26, axis = 22, ph = o.h - axis;
  const X = (ppm: number) => ((o.x0 - ppm) / (o.x0 - o.x1)) * o.w;
  const n = Math.max(200, Math.round(o.w * 2));
  const ppmPerPx = (o.x0 - o.x1) / o.w;
  const minW = ppmPerPx * 1.6;
  const lines = spectrumLines(sp, o.mhz);
  const ys = sampleSpectrum(lines, o.x0, o.x1, n, minW);
  const max = Math.max(...ys, 1e-9);
  const Y = (v: number) => ph - 2 - (v / max) * (ph - top);
  const path = (vals: number[]) => vals.map((v, i) => `${i ? 'L' : 'M'}${((i / (n - 1)) * o.w).toFixed(1)},${Y(v).toFixed(1)}`).join('');
  let out = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${o.w} ${o.h}" width="${o.w}" height="${o.h}" class="nmr-svg" font-family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Inter, sans-serif">`;
  if (c.bg) out += `<rect width="${o.w}" height="${o.h}" fill="${c.bg}"/>`;
  out += `<defs><clipPath id="${clip}"><rect x="0" y="0" width="${o.w}" height="${ph}"/></clipPath></defs>`;
  // Reference markers below the trace.
  for (const m of o.markers ?? []) {
    const col = m.kind === 'lit' ? c.lit : m.kind === 'solvent' ? c.muted : c.ref;
    if (m.lo !== undefined && m.hi !== undefined) {
      const a = X(Math.max(m.lo, m.hi)), b = X(Math.min(m.lo, m.hi));
      if (b < 0 || a > o.w) continue;
      out += `<rect x="${Math.max(0, a).toFixed(1)}" y="${top - 12}" width="${Math.max(2, Math.min(o.w, b) - Math.max(0, a)).toFixed(1)}" height="${ph - top + 12}" fill="${col}" opacity="0.09"/>`;
      out += `<rect x="${Math.max(0, a).toFixed(1)}" y="${ph - 5}" width="${Math.max(2, Math.min(o.w, b) - Math.max(0, a)).toFixed(1)}" height="4" rx="2" fill="${col}" opacity="0.6"/>`;
      if (markerText) {
        const mid = (Math.max(0, a) + Math.min(o.w, b)) / 2, tw = m.text.length * (fs - 1) * 0.55;
        out += `<text x="${mid.toFixed(1)}" y="${textY(mid - tw / 2, mid + tw / 2, top - 15)}" font-size="${fs - 1}" text-anchor="middle" fill="${col}">${esc(m.text)}</text>`;
      }
    } else if (m.v !== undefined) {
      const x = X(m.v);
      if (x < 0 || x > o.w) continue;
      out += `<line x1="${x.toFixed(1)}" x2="${x.toFixed(1)}" y1="${top - 10}" y2="${ph}" stroke="${col}" stroke-dasharray="3 3" stroke-width="1"/>`;
      if (markerText) {
        const tw = m.text.length * (fs - 1) * 0.55;
        const left = x + 3 + tw > o.w; // keep labels near the right edge inside the plot
        const a = left ? x - 3 - tw : x + 3;
        out += `<text x="${(left ? x - 3 : x + 3).toFixed(1)}" y="${textY(a, a + tw, top - 12)}" font-size="${fs - 1}" fill="${col}"${left ? ' text-anchor="end"' : ''}>${esc(m.text)}</text>`;
      }
    }
  }
  out += `<g clip-path="url(#${clip})">`;
  if (o.hl !== null && o.hl !== undefined) {
    const mine = sampleSpectrum(lines.filter((l) => l.peak === o.hl), o.x0, o.x1, n, minW);
    const col = PEAK_COLORS[o.hl % PEAK_COLORS.length];
    out += `<path d="${path(mine)}L${o.w},${ph - 2}L0,${ph - 2}Z" fill="${col}" opacity="0.18"/>`;
  }
  out += `<path d="${path(ys)}" fill="none" stroke="${c.fg}" stroke-width="1.2" stroke-linejoin="round"/></g>`;
  // Axis
  out += `<line x1="0" x2="${o.w}" y1="${ph}" y2="${ph}" stroke="${c.grid}"/>`;
  const step = tickStep(o.x0 - o.x1);
  const dec = step < 0.1 ? 2 : step < 1 ? 1 : 0;
  for (let t = Math.ceil(o.x1 / step) * step; t <= o.x0 + 1e-9; t += step) {
    const x = X(t);
    out += `<line x1="${x.toFixed(1)}" x2="${x.toFixed(1)}" y1="${ph}" y2="${ph + 4}" stroke="${c.muted}"/>`;
    if (x > 8 && x < o.w - 30) out += `<text x="${x.toFixed(1)}" y="${ph + 15}" font-size="${fs}" text-anchor="middle" fill="${c.muted}">${t.toFixed(dec)}</text>`;
  }
  out += `<text x="${o.w}" y="${o.h - 1}" font-size="${fs}" text-anchor="end" fill="${c.muted}">δ ppm</text>`;
  // Peak labels (staggered when crowded) and hover targets.
  const placed: { x: number; y: number }[] = [];
  const labels = sp.peaks.map((p, k) => {
    const ls = lines.filter((l) => l.peak === k);
    const lo = Math.min(...ls.map((l) => l.ppm)) - Math.max(p.width, minW), hi = Math.max(...ls.map((l) => l.ppm)) + Math.max(p.width, minW);
    const i0 = Math.max(0, Math.floor(((o.x0 - hi) / (o.x0 - o.x1)) * (n - 1))), i1 = Math.min(n - 1, Math.ceil(((o.x0 - lo) / (o.x0 - o.x1)) * (n - 1)));
    let m = 0; for (let i = i0; i <= i1; i++) m = Math.max(m, ys[i] ?? 0);
    return { k, x: X(p.shift), y: Y(m), a: X(hi), b: X(lo) };
  }).filter((l) => l.x >= -2 && l.x <= o.w + 2).sort((a, b) => a.x - b.x);
  for (const l of labels) {
    let y = Math.max(top - 2, l.y - 10);
    while (placed.some((p) => Math.abs(p.x - l.x) < 15 && Math.abs(p.y - y) < 14)) y -= 14;
    y = Math.max(8, y);
    placed.push({ x: l.x, y });
    const col = PEAK_COLORS[l.k % PEAK_COLORS.length];
    const p = sp.peaks[l.k];
    out += `<g class="nmr-pk" data-peak="${l.k}"><rect x="${(Math.min(l.a, l.x - 6)).toFixed(1)}" y="0" width="${Math.max(12, l.b - l.a).toFixed(1)}" height="${ph}" fill="transparent"/>`;
    out += `<circle cx="${l.x.toFixed(1)}" cy="${y.toFixed(1)}" r="7" fill="${col}"/><text x="${l.x.toFixed(1)}" y="${(y + 3.4).toFixed(1)}" font-size="${p.label.length > 1 ? 8 : 9.5}" font-weight="600" text-anchor="middle" fill="#fff">${p.label}</text></g>`;
  }
  return out + '</svg>';
}

function markersFor(n: Nucleus): Marker[] {
  const out: Marker[] = [];
  const pick = (vs: RefValue[] | undefined, text: string, kind: Marker['kind']) => {
    for (const v of vs ?? []) out.push(Array.isArray(v) ? { lo: v[0], hi: v[1], text, kind } : { v, text: `${text} ${v.toFixed(n === '1H' ? 2 : 1)}`, kind });
  };
  if (st.showSolvent) {
    const s = SOLVENT_PEAKS[st.solvent];
    pick(n === '1H' ? s.h : s.c, s.label, 'solvent');
    if (n === '1H' && st.solvent !== 'D2O') pick([s.water], 'H₂O', 'solvent');
  }
  if (st.showLit && match) pick(n === '1H' ? match.h : match.c, 'lit.', 'lit');
  for (const id of st.pinned) {
    const r = NMR_REFS.find((x) => x.id === id);
    if (r) pick(n === '1H' ? r.h : r.c, shortName(r), 'ref');
  }
  return out;
}

const shortName = (r: NmrRef) => r.name.replace(/\s*\(.*\)$/, '').replace(/ in .*$/, '');

// ---------- inspector panel ----------

function controls(): string {
  return `<div class="nmr-controls">
    <div class="seg" role="group" aria-label="Nucleus">${(['1H', '13C'] as Nucleus[]).map((n) => `<button data-nuc="${n}" class="${st.nucleus === n ? 'on' : ''}" aria-pressed="${st.nucleus === n}">${NUC_LABEL[n]}</button>`).join('')}</div>
    <select data-nmr="solvent" title="Solvent">${NMR_SOLVENTS.map((s) => `<option value="${s}" ${s === st.solvent ? 'selected' : ''}>${SOLVENT_PEAKS[s].label}</option>`).join('')}</select>
    <select data-nmr="mhz" title="Spectrometer frequency (¹H)">${[300, 400, 500, 600, 800].map((m) => `<option value="${m}" ${m === st.mhz ? 'selected' : ''}>${m} MHz</option>`).join('')}</select>
  </div>`;
}

function wireControls(root: HTMLElement, rerender: () => void) {
  for (const b of root.querySelectorAll<HTMLButtonElement>('[data-nuc]')) b.addEventListener('click', () => { st.nucleus = b.dataset.nuc as Nucleus; persist(); rerender(); });
  root.querySelector<HTMLSelectElement>('[data-nmr="solvent"]')?.addEventListener('change', (e) => {
    st.solvent = (e.target as HTMLSelectElement).value as NmrSolvent; persist();
    for (const k of Object.keys(results)) delete results[k as Nucleus];
    renderPanel(); if (dialog?.open) renderDialog();
  });
  root.querySelector<HTMLSelectElement>('[data-nmr="mhz"]')?.addEventListener('change', (e) => {
    st.mhz = Number((e.target as HTMLSelectElement).value); persist(); renderPanel(); if (dialog?.open) renderDialog();
  });
}

function peakTable(sp: NmrSpectrum, compact: boolean): string {
  const n = sp.nucleus;
  const head = n === '1H' ? '<th></th><th>δ (ppm)</th><th>mult.</th><th>∫</th><th>assignment</th>' : '<th></th><th>δ (ppm)</th><th>C</th><th>assignment</th>';
  const rows = sp.peaks.map((p, k) => {
    const col = PEAK_COLORS[k % PEAK_COLORS.length];
    const dot = `<span class="pk-dot" style="background:${col}">${p.label}</span>`;
    const shift = `${p.rough ? '~' : ''}${formatShift(p, n)}`;
    const j = !compact && n === '1H' && p.couplings.length ? `<span class="muted"> ${couplingText(p)}</span>` : '';
    const what = `${esc(p.assignment)}${p.exchangeable ? ' <span class="muted">(exch.)</span>' : ''}${n === '13C' && p.quaternary ? ' <span class="muted">(q)</span>' : ''}`;
    return n === '1H'
      ? `<tr data-peak="${k}"><td>${dot}</td><td>${shift}</td><td>${p.mult}${j}</td><td>${p.count}H</td><td>${what}</td></tr>`
      : `<tr data-peak="${k}"><td>${dot}</td><td>${shift}</td><td>${p.count > 1 ? `×${p.count}` : ''}</td><td>${what}</td></tr>`;
  }).join('');
  return `<table class="nmr-table"><thead><tr>${head}</tr></thead><tbody>${rows}</tbody></table>`;
}

function method(sp: NmrSpectrum): string {
  const n = sp.nucleus;
  const base = n === '1H'
    ? 'Additive increments (Curphy–Morrison, Pascual–Meier–Simon, Pretsch), first-order couplings. Typically ±0.2–0.3 ppm.'
    : 'Lindeman–Adams skeleton plus Pretsch substituent increments, ¹H-decoupled. Typically ±3–5 ppm.';
  const poly = sp.isPolymer ? ' Per repeat unit, as part of an endless chain (the two * ends joined); polymer signals are broad and tacticity splitting is not predicted.' : '';
  return `<div class="muted method">${base}${poly} ${sp.notes.map(esc).join(' ')}</div>`;
}

function litLine(n: Nucleus): string {
  if (!match) return '';
  const vs = n === '1H' ? match.h : match.c;
  if (!vs?.length) return '';
  return `<div class="nmr-lit small"><span class="lit-dot"></span>Literature, ${esc(match.name)}: ${fmtRef(vs, n === '1H' ? 2 : 1)} <span class="muted">(${esc(match.note ?? '')})</span></div>`;
}

let hlPanel: number | null = null;
export function renderPanel() {
  const el = document.getElementById('nmr');
  if (!el || !st.open) return;
  const btn = '<button class="mini wide" id="open-nmr">Open spectrum &amp; reference peaks…</button>';
  if (!doc) {
    el.innerHTML = `<div class="muted small">Draw a molecule or a polymer repeat unit to see its predicted ¹H and ¹³C NMR spectrum.</div>${btn}`;
  } else {
    const r = result(st.nucleus);
    if (!r.ok) el.innerHTML = `${controls()}<div class="muted small">Not available: ${esc(r.reason)}</div>${btn}`;
    else {
      const sp = r.value;
      const [x0, x1] = autoRange(sp);
      el.innerHTML = `${controls()}<div class="nmr-plot">${plotSvg(sp, { w: 256, h: 132, x0, x1, mhz: st.mhz, hl: hlPanel, markers: markersFor(st.nucleus).filter((m) => m.kind === 'solvent'), markerText: false })}</div>
        ${sp.peaks.length ? peakTable(sp, true) : `<div class="muted small">No ${NUC_LABEL[sp.nucleus]} signals.</div>`}
        ${litLine(sp.nucleus)}${method(sp)}${btn}`;
      wireHover(el, sp, (k) => { hlPanel = k; renderPanelPlot(); });
    }
    wireControls(el, renderPanel);
  }
  el.querySelector('#open-nmr')?.addEventListener('click', openNmrDialog);
}

function renderPanelPlot() {
  const el = document.getElementById('nmr');
  const plot = el?.querySelector('.nmr-plot');
  const r = result(st.nucleus);
  if (!plot || !r.ok) return;
  const [x0, x1] = autoRange(r.value);
  plot.innerHTML = plotSvg(r.value, { w: 256, h: 132, x0, x1, mhz: st.mhz, hl: hlPanel, markers: markersFor(st.nucleus).filter((m) => m.kind === 'solvent'), markerText: false });
  for (const row of el!.querySelectorAll<HTMLElement>('tr[data-peak]')) row.classList.toggle('hl', Number(row.dataset.peak) === hlPanel);
}

/** Hovering a peak (in the plot or the table) rings its atoms on the canvas. */
function wireHover(root: HTMLElement, sp: NmrSpectrum, onPeak: (k: number | null) => void) {
  let cur: number | null = null;
  const set = (k: number | null) => {
    if (k === cur) return;
    cur = k;
    onPeak(k);
    if (k === null) { host.setMarks(null); return; }
    const col = PEAK_COLORS[k % PEAK_COLORS.length];
    host.setMarks(new Map(sp.peaks[k].atoms.map((id) => [id, col])));
  };
  root.addEventListener('pointerover', (e) => {
    const t = (e.target as Element).closest<HTMLElement | SVGElement>('[data-peak]');
    if (t && root.contains(t)) set(Number((t as HTMLElement).dataset.peak));
  });
  root.addEventListener('pointerleave', () => set(null));
  root.addEventListener('pointerout', (e) => {
    const to = (e as PointerEvent).relatedTarget as Element | null;
    if (!to || !to.closest('[data-peak]')) set(null);
  });
}

// ---------- dialog ----------

export function openNmrDialog() {
  if (!dialog) {
    dialog = document.createElement('dialog');
    dialog.className = 'analyzer nmr-dialog';
    dialog.setAttribute('aria-label', 'NMR spectrum');
    document.body.appendChild(dialog);
    dialog.addEventListener('click', (e) => { if (e.target === dialog) dialog!.close(); });
    dialog.addEventListener('close', () => host.setMarks(null));
  }
  renderDialog();
  dialog.showModal();
}

const hlDlg: Partial<Record<Nucleus, number | null>> = {};

function dialogPlot(n: Nucleus, w: number, h: number): string {
  const r = result(n);
  if (!r.ok) return `<p class="muted small">${NUC_LABEL[n]}: not available. ${esc(r.reason)}</p>`;
  const [x0, x1] = zoom[n] ?? autoRange(r.value);
  return plotSvg(r.value, { w, h, x0, x1, mhz: st.mhz, hl: hlDlg[n] ?? null, markers: markersFor(n), font: 11 });
}

function structureMap(): string {
  if (!doc) return '';
  const r = result(st.nucleus);
  if (!r.ok) return '';
  const byAtom = new Map<number, { t: string; c: string }[]>();
  r.value.peaks.forEach((p, k) => { for (const id of p.atoms) byAtom.set(id, [...(byAtom.get(id) ?? []), { t: p.label, c: PEAK_COLORS[k % PEAK_COLORS.length] }]); });
  let svg = standaloneSvg(doc, { theme: host.isDark() ? DARK : LIGHT, colorAtoms: true, pad: 22 });
  // Put each atom's peak letters just outside the structure: opposite the mean direction of its
  // bonds (so they never sit on a bond or inside a ring), centred on that point.
  const pos = new Map(doc.atoms.map((a) => [a.id, a]));
  const halo = host.isDark() ? '#1f1f22' : '#fff';
  let tags = '';
  for (const a of doc.atoms) {
    const ls = byAtom.get(a.id);
    if (!ls) continue;
    let dx = 0, dy = 0;
    for (const b of doc.bonds) {
      const o = b.a === a.id ? pos.get(b.b) : b.b === a.id ? pos.get(b.a) : undefined;
      if (!o) continue;
      const len = Math.hypot(o.x - a.x, o.y - a.y) || 1;
      dx += (o.x - a.x) / len; dy += (o.y - a.y) / len;
    }
    let ux = -dx, uy = -dy;
    const m = Math.hypot(ux, uy);
    if (m < 0.3) {
      // Linear or isolated atom: go perpendicular to the bonds (or straight up).
      const n = Math.hypot(dx, dy);
      if (n > 0.01) { ux = -dy / n; uy = dx / n; } else { ux = 0; uy = -1; }
      if (uy > 0) { ux = -ux; uy = -uy; }
    } else { ux /= m; uy /= m; }
    // A heteroatom's own label (OH, NH2, HO…) runs sideways, so put its letters above or below it.
    if (a.el !== 'C') { uy = uy > 0 ? 1 : -1; ux = 0; }
    const text = ls.map((l) => l.t).join(',');
    const r = (a.el === 'C' ? 12 : 16) + Math.abs(ux) * Math.max(0, text.length - 1) * 3.5;
    const spans = ls.map((l, i) => `${i ? `<tspan fill="${halo === '#fff' ? '#888' : '#aaa'}">,</tspan>` : ''}<tspan fill="${l.c}">${l.t}</tspan>`).join('');
    tags += `<text x="${(a.x + ux * r).toFixed(1)}" y="${(a.y + uy * r).toFixed(1)}" text-anchor="middle" dominant-baseline="central" font-size="13" font-weight="700" stroke="${halo}" stroke-width="3" paint-order="stroke" font-family="-apple-system, sans-serif">${spans}</text>`;
  }
  svg = svg.replace(/<\/svg>$/, `${tags}</svg>`);
  return svg;
}

function refList(): string {
  const q = st.query.trim().toLowerCase();
  const items = NMR_REFS.filter((r) => (st.cat === 'all' || r.cat === st.cat) && (!q || `${r.name} ${r.note ?? ''} ${fmtRef(r.h, 2)} ${fmtRef(r.c, 1)}`.toLowerCase().includes(q)));
  if (match && !items.includes(match) && (st.cat === 'all' || st.cat === match.cat) && !q) items.unshift(match);
  if (!items.length) return '<p class="muted small">No reference matches. Try a group (ester, aromatic), a solvent or a polymer name, or a shift such as 3.6.</p>';
  return items.map((r) => {
    const on = st.pinned.includes(r.id);
    const isMatch = r === match;
    return `<button class="ref-item ${on ? 'on' : ''}" data-ref="${r.id}" aria-pressed="${on}">
      <span class="ref-name">${esc(r.name)}${isMatch ? ' <span class="ref-badge">matches drawing</span>' : ''}</span>
      ${r.h?.length ? `<span class="ref-vals"><b>¹H</b> ${fmtRef(r.h, 2)}</span>` : ''}
      ${r.c?.length ? `<span class="ref-vals"><b>¹³C</b> ${fmtRef(r.c, 1)}</span>` : ''}
      ${r.note ? `<span class="muted ref-note">${esc(r.note)}</span>` : ''}
    </button>`;
  }).join('');
}

function renderDialog() {
  const d = dialog!;
  const refScroll = d.querySelector('.nmr-refs')?.scrollTop ?? 0;
  const hadFocus = document.activeElement?.id === 'nmr-search';
  const title = doc ? (label ? esc(label) : 'drawn structure') : 'no structure';
  const anyOk = (['1H', '13C'] as Nucleus[]).some((n) => result(n).ok);
  const polymer = (['1H', '13C'] as Nucleus[]).map((n) => result(n)).find((r) => r.ok && r.value.isPolymer);
  d.innerHTML = `
    <header class="an-head">
      <h2>Predicted NMR <span class="muted nmr-sub">${title}${polymer ? ' · per repeat unit' : ''}</span></h2>
      <div class="nmr-head-actions">
        <button class="mini" data-act="csv" ${anyOk ? '' : 'disabled'}>Export CSV</button>
        <button class="mini" data-act="png" ${anyOk ? '' : 'disabled'}>Export PNG</button>
        <button class="mini" data-act="close" aria-label="Close">Close</button>
      </div>
    </header>
    <div class="nmr-body">
      <div class="nmr-main">
        <div class="nmr-toolbar">${controls()}
          <button class="chip ${st.showSolvent ? 'on' : ''}" data-chip="solvent" aria-pressed="${st.showSolvent}">Solvent &amp; water</button>
          ${match ? `<button class="chip ${st.showLit ? 'on' : ''}" data-chip="lit" aria-pressed="${st.showLit}">Literature</button>` : ''}
          <span class="muted small">Drag across a spectrum to zoom, double-click to reset.</span>
        </div>
        ${doc ? '' : '<p class="muted small">Draw a structure on the canvas to predict its spectrum. The reference list works without one.</p>'}
        ${doc ? `<div class="nmr-plot-wrap" data-plot="1H"><div class="nmr-plot-title">¹H <span class="muted">${st.mhz} MHz, ${SOLVENT_PEAKS[st.solvent].label}</span></div>${dialogPlot('1H', 760, 250)}</div>
        <div class="nmr-plot-wrap" data-plot="13C"><div class="nmr-plot-title">¹³C{¹H} <span class="muted">${Math.round(st.mhz * 0.2515)} MHz</span></div>${dialogPlot('13C', 760, 170)}</div>
        <div class="nmr-tables">
          <div class="nmr-structure">${structureMap()}</div>
          ${(['1H', '13C'] as Nucleus[]).map((n) => { const r = result(n); return r.ok && r.value.peaks.length ? `<div data-table="${n}"><div class="nmr-plot-title">${NUC_LABEL[n]} peaks</div>${peakTable(r.value, false)}</div>` : ''; }).join('')}
        </div>
        ${(() => { const r = result(st.nucleus); return r.ok ? method(r.value) : ''; })()}` : ''}
      </div>
      <aside class="nmr-side">
        <h3>Reference peaks</h3>
        <input id="nmr-search" class="nmr-search" type="search" placeholder="Search: ester, aromatic, DMSO, PMMA, 3.6…" value="${esc(st.query)}" />
        <div class="nmr-cats">${REF_CATEGORIES.map((c) => `<button class="chip ${st.cat === c.id ? 'on' : ''}" data-cat="${c.id}">${c.label}</button>`).join('')}</div>
        <div class="nmr-refs">${refList()}</div>
        <p class="muted small an-note">Click a reference to pin it on the spectra: single values as dashed lines, ranges as shaded bands. Hover a peak to see its atoms.</p>
      </aside>
    </div>`;
  d.querySelector('.nmr-refs')!.scrollTop = refScroll;
  if (hadFocus) { const s = d.querySelector<HTMLInputElement>('#nmr-search')!; s.focus(); s.setSelectionRange(s.value.length, s.value.length); }
  wireDialog(d);
}

function rerenderPlots() {
  const d = dialog!;
  for (const n of ['1H', '13C'] as Nucleus[]) {
    const wrap = d.querySelector<HTMLElement>(`[data-plot="${n}"]`);
    const svg = wrap?.querySelector('svg');
    if (wrap && svg) svg.outerHTML = dialogPlot(n, 760, n === '1H' ? 250 : 170);
    for (const row of d.querySelectorAll<HTMLElement>(`[data-table="${n}"] tr[data-peak]`)) row.classList.toggle('hl', Number(row.dataset.peak) === hlDlg[n]);
  }
  wireZoom(d);
}

function wireDialog(d: HTMLDialogElement) {
  d.querySelector('[data-act="close"]')?.addEventListener('click', () => d.close());
  d.querySelector('[data-act="csv"]')?.addEventListener('click', exportCsv);
  d.querySelector('[data-act="png"]')?.addEventListener('click', () => { exportPng().catch(() => host.notify('Could not export the spectrum image.')); });
  wireControls(d, () => { renderPanel(); renderDialog(); });
  for (const b of d.querySelectorAll<HTMLButtonElement>('[data-chip]')) b.addEventListener('click', () => {
    if (b.dataset.chip === 'solvent') st.showSolvent = !st.showSolvent; else st.showLit = !st.showLit;
    persist(); renderDialog(); renderPanel();
  });
  for (const b of d.querySelectorAll<HTMLButtonElement>('[data-cat]')) b.addEventListener('click', () => { st.cat = b.dataset.cat as State['cat']; persist(); renderDialog(); });
  d.querySelector<HTMLInputElement>('#nmr-search')?.addEventListener('input', (e) => {
    st.query = (e.target as HTMLInputElement).value; persist();
    d.querySelector('.nmr-refs')!.innerHTML = refList(); wireRefs(d);
  });
  wireRefs(d);
  for (const n of ['1H', '13C'] as Nucleus[]) {
    const r = result(n);
    if (!r.ok) continue;
    const parts = [d.querySelector<HTMLElement>(`[data-plot="${n}"]`), d.querySelector<HTMLElement>(`[data-table="${n}"]`)].filter((x): x is HTMLElement => !!x);
    for (const p of parts) wireHover(p, r.value, (k) => { hlDlg[n] = k; rerenderPlots(); });
  }
  wireZoom(d);
}

function wireRefs(d: HTMLDialogElement) {
  for (const b of d.querySelectorAll<HTMLButtonElement>('[data-ref]')) b.addEventListener('click', () => {
    const id = b.dataset.ref!;
    st.pinned = st.pinned.includes(id) ? st.pinned.filter((x) => x !== id) : [...st.pinned, id];
    persist();
    b.classList.toggle('on'); b.setAttribute('aria-pressed', String(st.pinned.includes(id)));
    rerenderPlots();
  });
}

/** Drag across a plot to zoom into that ppm window; double-click to reset. */
function wireZoom(d: HTMLDialogElement) {
  for (const n of ['1H', '13C'] as Nucleus[]) {
    const wrap = d.querySelector<HTMLElement>(`[data-plot="${n}"]`);
    const svg = wrap?.querySelector<SVGSVGElement>('svg');
    const r = result(n);
    if (!wrap || !svg || !r.ok || svg.dataset.zoomWired) continue;
    svg.dataset.zoomWired = '1';
    const range = () => zoom[n] ?? autoRange(r.value);
    const toPpm = (clientX: number) => { const b = svg.getBoundingClientRect(); const [x0, x1] = range(); return x0 - ((clientX - b.left) / b.width) * (x0 - x1); };
    let start: number | null = null;
    let box: HTMLDivElement | null = null;
    svg.addEventListener('pointerdown', (e) => {
      start = e.clientX;
      svg.setPointerCapture(e.pointerId);
      box = document.createElement('div'); box.className = 'nmr-zoombox'; wrap.appendChild(box);
    });
    svg.addEventListener('pointermove', (e) => {
      if (start === null || !box) return;
      const b = svg.getBoundingClientRect(), w = wrap.getBoundingClientRect();
      box.style.left = `${Math.min(start, e.clientX) - w.left}px`; box.style.width = `${Math.abs(e.clientX - start)}px`;
      box.style.top = `${b.top - w.top}px`; box.style.height = `${b.height}px`;
    });
    svg.addEventListener('pointerup', (e) => {
      box?.remove(); box = null;
      if (start === null) return;
      const a = toPpm(start), b = toPpm(e.clientX);
      start = null;
      const span = Math.abs(a - b);
      if (span < (n === '1H' ? 0.02 : 0.5)) return;
      zoom[n] = [Math.max(a, b), Math.min(a, b)];
      rerenderPlots();
    });
    svg.addEventListener('dblclick', () => { delete zoom[n]; rerenderPlots(); });
  }
}

// ---------- export ----------

function exportCsv() {
  const rows = ['nucleus,label,shift_ppm,range_ppm,multiplicity,J_Hz,count,assignment,estimate_quality'];
  for (const n of ['1H', '13C'] as Nucleus[]) {
    const r = result(n);
    if (!r.ok) continue;
    for (const p of r.value.peaks) {
      rows.push([n, p.label, p.shift.toFixed(n === '1H' ? 2 : 1), p.range ? `${p.range[0]}-${p.range[1]}` : '', p.mult,
        p.couplings.map((c) => c.j.toFixed(1)).join(' '), p.count, `"${p.assignment.replace(/"/g, '""')}"`, p.rough ? 'rough' : 'normal'].join(','));
    }
  }
  saveFile('predicted-nmr.csv', rows.join('\n') + '\n', { name: 'CSV', extensions: ['csv'] }).catch(() => host.notify('Could not save the file.'));
}

async function exportPng() {
  const W = 1000;
  const colors: Colors = { fg: '#1d1d1f', muted: '#6e6e73', grid: '#d2d2d7', ref: '#8e44ad', lit: '#1e9e3a', band: '#8e44ad', bg: '#ffffff' };
  const parts: string[] = [];
  let y = 0;
  const head = (t: string) => { parts.push(`<text x="0" y="${y + 16}" font-size="15" font-weight="600" fill="#1d1d1f" font-family="-apple-system, sans-serif">${esc(t)}</text>`); y += 26; };
  head(`Predicted NMR · ${label || 'structure'} (MolCanvas estimate)`);
  for (const n of ['1H', '13C'] as Nucleus[]) {
    const r = result(n);
    if (!r.ok) continue;
    const [x0, x1] = zoom[n] ?? autoRange(r.value);
    const h = n === '1H' ? 300 : 220;
    head(n === '1H' ? `¹H, ${st.mhz} MHz, ${SOLVENT_PEAKS[st.solvent].label}` : `¹³C{¹H}, ${Math.round(st.mhz * 0.2515)} MHz`);
    parts.push(`<g transform="translate(0 ${y})">${plotSvg(r.value, { w: W, h, x0, x1, mhz: st.mhz, markers: markersFor(n), colors, font: 12 })}</g>`);
    y += h + 16;
  }
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${W + 40}" height="${y + 20}" viewBox="-20 -10 ${W + 40} ${y + 20}"><rect x="-20" y="-10" width="${W + 40}" height="${y + 20}" fill="#fff"/>${parts.join('')}</svg>`;
  const img = new Image();
  img.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}`;
  await img.decode();
  const scale = 2;
  const canvas = document.createElement('canvas');
  canvas.width = (W + 40) * scale; canvas.height = (y + 20) * scale;
  const ctx = canvas.getContext('2d')!;
  ctx.scale(scale, scale);
  ctx.drawImage(img, 0, 0);
  const blob = await new Promise<Blob | null>((res) => canvas.toBlob(res, 'image/png'));
  if (!blob) throw new Error('no image');
  await saveFile('predicted-nmr.png', new Uint8Array(await blob.arrayBuffer()), { name: 'PNG image', extensions: ['png'] });
}

