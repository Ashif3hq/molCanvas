import './style.css';
import { Editor, Tool } from './editor';
import { fromSmiles, isElement, properties, toMolfile, toSmiles } from './chem';
import { ICONS } from './icons';
import { CATEGORIES, Preset, PRESETS } from './presets';
import { DARK, LIGHT, standaloneSvg, atomColor } from './render';
import { Doc, emptyDoc } from './model';
import { orientPolymer } from './geometry';
import { isDesktop, isMac, openTextFile, saveFile } from './platform';

const $ = <T extends HTMLElement = HTMLElement>(id: string) => document.getElementById(id) as T;
const mod = isMac ? '⌘' : 'Ctrl+';
const STORAGE_KEY = 'molcanvas.doc.v1';

// ---------- theme ----------

const prefersDark = window.matchMedia('(prefers-color-scheme: dark)');
let themeChoice: 'light' | 'dark' | null = null;
try { themeChoice = localStorage.getItem('molcanvas.theme') as 'light' | 'dark' | null; } catch { /* storage unavailable */ }
const isDark = () => (themeChoice ?? (prefersDark.matches ? 'dark' : 'light')) === 'dark';
function applyTheme() {
  document.documentElement.dataset.theme = isDark() ? 'dark' : 'light';
  if (editor) { editor.theme = isDark() ? DARK : LIGHT; editor.render(); }
  const btn = document.getElementById('theme-btn');
  if (btn) btn.innerHTML = isDark() ? ICONS.sun : ICONS.moon;
}
prefersDark.addEventListener('change', applyTheme);

if (isDesktop() && isMac) document.documentElement.classList.add('desktop-mac');

// ---------- editor ----------

let editor: Editor;
applyTheme();
editor = new Editor($('canvas') as unknown as SVGSVGElement, $('canvas-wrap'), isDark() ? DARK : LIGHT);
editor.isElement = isElement;

// Group abbreviations for the text tool (Me, OMe, CO2H, Ph, Boc …)
const ABBREV = new Map<string, Preset>();
for (const p of PRESETS) {
  if (p.kind !== 'group') continue;
  for (const tok of [p.name, ...(p.aka ?? '').split(/\s+/)]) {
    if (tok && /[A-Z]/.test(tok) && tok.length <= 8 && !ABBREV.has(tok)) ABBREV.set(tok, p);
  }
}
const fragCache = new Map<string, Doc>();
function presetFrag(p: Preset): Doc {
  let d = fragCache.get(p.smiles);
  if (!d) {
    d = fromSmiles(p.smiles);
    if (p.kind === 'polymer') orientPolymer(d);
    fragCache.set(p.smiles, d);
  }
  return d;
}
editor.abbreviation = (label) => {
  const p = ABBREV.get(label);
  return p ? presetFrag(p) : null;
};

// ---------- toast & status ----------

let toastTimer = 0;
function toast(msg: string) {
  const el = $('toast');
  el.textContent = msg;
  el.classList.add('show');
  clearTimeout(toastTimer);
  toastTimer = window.setTimeout(() => el.classList.remove('show'), 2600);
}
// ---------- info bar: explains whatever is under the pointer ----------

let toolHint = { title: '', detail: '' };
let warnUntil = 0;
function showInfo(title: string, detail: string, warn = false) {
  $('info-title').textContent = title;
  $('info-detail').textContent = detail;
  $('status').classList.toggle('warn', warn);
}
function splitHint(msg: string) {
  const i = msg.indexOf(':');
  return i > 0 && i < 40 ? { title: msg.slice(0, i), detail: msg.slice(i + 1).trim() } : { title: 'Tip', detail: msg };
}
editor.onStatus = (msg) => {
  toolHint = splitHint(msg);
  if (Date.now() > warnUntil) showInfo(toolHint.title, toolHint.detail);
};
editor.onHover = (info) => {
  if (Date.now() < warnUntil && !info?.warn) return;
  if (info) showInfo(info.title, info.detail, !!info.warn);
  else showInfo(toolHint.title, toolHint.detail);
};
editor.onWarn = (msg) => {
  warnUntil = Date.now() + 3500;
  showInfo('Not allowed', msg, true);
  toast(msg);
};
/** Show a description while the pointer is over a toolbar button or library card. */
function explainOnHover(el: HTMLElement, title: string, detail: string) {
  el.addEventListener('mouseenter', () => showInfo(title, detail));
  el.addEventListener('mouseleave', () => showInfo(toolHint.title, toolHint.detail));
}

// ---------- top bar ----------

function iconButton(icon: string, title: string, onClick: (e: MouseEvent) => void, label?: string): HTMLButtonElement {
  const b = document.createElement('button');
  b.className = 'icon-btn';
  b.innerHTML = ICONS[icon] + (label ? `<span class="lbl">${label}</span>` : '');
  b.title = title;
  b.addEventListener('click', onClick);
  return b;
}

const fileActions = $('file-actions');
fileActions.append(
  iconButton('library', 'Show / hide library', () => $('app').classList.toggle('lib-hidden')),
  iconButton('new', `New drawing (${mod}N)`, () => newDrawing()),
  iconButton('open', `Open .mcv, .mol, .sdf or .smi (${mod}O)`, () => openFile()),
  iconButton('save', `Save as .mcv (${mod}S)`, () => saveDrawing()),
  iconButton('export', 'Export PNG, SVG, MOL or SMILES', (e) => showExportMenu(e.currentTarget as HTMLElement), 'Export'),
  iconButton('copy', `Copy as image for Word, PowerPoint, Keynote (${mod}⇧C)`, () => copyImage()),
);

const undoBtn = iconButton('undo', `Undo (${mod}Z)`, () => editor.undo());
const redoBtn = iconButton('redo', `Redo (${mod}⇧Z)`, () => editor.redo());
$('edit-actions').append(undoBtn, redoBtn);

$('transform-actions').append(
  iconButton('clean', 'Clean up structure layout', () => editor.cleanUp(), 'Clean'),
  iconButton('rotate', 'Rotate 15° (selection or all)', () => editor.rotate(15)),
  iconButton('flipH', 'Flip horizontal', () => editor.flip('h')),
  iconButton('flipV', 'Flip vertical', () => editor.flip('v')),
);

const colorBtn = iconButton('palette', 'Colour heteroatoms', () => {
  editor.colorAtoms = !editor.colorAtoms;
  colorBtn.classList.toggle('active', editor.colorAtoms);
  editor.render();
  renderLibrary();
});
colorBtn.classList.add('active');
const themeBtn = iconButton('moon', 'Toggle dark mode', () => {
  themeChoice = isDark() ? 'light' : 'dark';
  try { localStorage.setItem('molcanvas.theme', themeChoice); } catch { /* ignore */ }
  applyTheme();
  renderLibrary();
});
themeBtn.id = 'theme-btn';
$('view-actions').append(
  iconButton('zoomOut', `Zoom out (${mod}-)`, () => editor.zoomBy(1 / 1.25)),
  iconButton('fit', `Fit to window (${mod}0)`, () => editor.fitView()),
  iconButton('zoomIn', `Zoom in (${mod}+)`, () => editor.zoomBy(1.25)),
  colorBtn,
  themeBtn,
);
applyTheme();

// ---------- tool rail ----------

interface ToolDef { key: string; icon?: string; label?: string; title: string; tool: Tool; cls?: string; }
const TOOL_DEFS: (ToolDef | 'sep')[] = [
  { key: 'select', icon: 'select', title: 'Select & move (V)', tool: { kind: 'select' } },
  { key: 'erase', icon: 'erase', title: 'Eraser (E)', tool: { kind: 'erase' } },
  'sep',
  { key: 'b1', icon: 'single', title: 'Single bond (1 / B). Click a bond to cycle single → double → triple', tool: { kind: 'bond', order: 1, stereo: 'none' } },
  { key: 'b2', icon: 'double', title: 'Double bond (2)', tool: { kind: 'bond', order: 2, stereo: 'none' } },
  { key: 'b3', icon: 'triple', title: 'Triple bond (3)', tool: { kind: 'bond', order: 3, stereo: 'none' } },
  { key: 'wedge', icon: 'wedge', title: 'Wedge bond (up). Click again to flip', tool: { kind: 'bond', order: 1, stereo: 'wedge' } },
  { key: 'hash', icon: 'hash', title: 'Hashed bond (down). Click again to flip', tool: { kind: 'bond', order: 1, stereo: 'hash' } },
  { key: 'wavy', icon: 'wavy', title: 'Wavy bond (unknown stereo)', tool: { kind: 'bond', order: 1, stereo: 'wavy' } },
  { key: 'chain', icon: 'chain', title: 'Carbon chain (G): drag to set length', tool: { kind: 'chain' } },
  'sep',
  { key: 'benzene', icon: 'benzene', title: 'Benzene ring (R)', tool: { kind: 'ring', size: 6, aromatic: true } },
  { key: 'r3', icon: 'ring3', title: 'Cyclopropane', tool: { kind: 'ring', size: 3, aromatic: false } },
  { key: 'r4', icon: 'ring4', title: 'Cyclobutane', tool: { kind: 'ring', size: 4, aromatic: false } },
  { key: 'r5', icon: 'ring5', title: 'Cyclopentane', tool: { kind: 'ring', size: 5, aromatic: false } },
  { key: 'r6', icon: 'ring6', title: 'Cyclohexane', tool: { kind: 'ring', size: 6, aromatic: false } },
  { key: 'r7', icon: 'ring7', title: 'Cycloheptane', tool: { kind: 'ring', size: 7, aromatic: false } },
  { key: 'r8', icon: 'ring8', title: 'Cyclooctane', tool: { kind: 'ring', size: 8, aromatic: false } },
  'sep',
  ...['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'].map((el): ToolDef => ({ key: `el-${el}`, label: el, cls: 'el', title: `${el} atom (hover an atom and press ${el === 'Cl' ? 'L' : el === 'Br' ? 'B' : el})`, tool: { kind: 'atom', el } })),
  { key: 'el-*', icon: 'polymer', title: 'Polymer end / attachment point (*): drawn as a repeat-unit bracket', tool: { kind: 'atom', el: '*' } },
  { key: 'ptable', icon: 'table', title: 'Periodic table', tool: { kind: 'atom', el: 'C' } },
  'sep',
  { key: 'plus', icon: 'plus', title: 'Positive charge (hover + "+")', tool: { kind: 'charge', delta: 1 } },
  { key: 'minus', icon: 'minus', title: 'Negative charge (hover + "−")', tool: { kind: 'charge', delta: -1 } },
  'sep',
  { key: 'arrow', icon: 'arrow', title: 'Reaction arrow (A)', tool: { kind: 'arrow', arrow: 'forward' } },
  { key: 'equilibrium', icon: 'equilibrium', title: 'Equilibrium arrow', tool: { kind: 'arrow', arrow: 'equilibrium' } },
  { key: 'retro', icon: 'retro', title: 'Retrosynthesis arrow', tool: { kind: 'arrow', arrow: 'retro' } },
  { key: 'resonance', icon: 'resonance', title: 'Resonance arrow', tool: { kind: 'arrow', arrow: 'resonance' } },
  { key: 'text', icon: 'text', title: 'Text & labels (T). Click an atom to type an element or group (OMe, CO2H, Boc…)', tool: { kind: 'text' } },
];

const TOOL_EXPLAIN: Record<string, string> = {
  select: 'Click to select, drag to box-select, drag a selection to move it, double-click to select a whole molecule. Shortcut V.',
  erase: 'Click or drag over atoms, bonds, arrows and text to delete them. Shortcut E.',
  b1: 'Single (σ) bond. Click an atom to grow a bond, drag to aim it, click an existing bond to cycle single → double → triple.',
  b2: 'Double bond: one σ and one π bond, as in alkenes (C=C) and carbonyls (C=O).',
  b3: 'Triple bond: one σ and two π bonds, as in alkynes (C≡C) and nitriles (C≡N). The atoms become linear (sp).',
  wedge: 'Stereo bond pointing towards you. The narrow end sits on the stereocentre. Click a wedge again to flip it.',
  hash: 'Stereo bond pointing away from you. The narrow end sits on the stereocentre. Click again to flip it.',
  wavy: 'Bond with unknown or mixed stereochemistry (e.g. a racemic centre or E/Z mixture).',
  chain: 'Drag to draw a zig-zag carbon chain; the bar shows how many bonds. Shortcut G.',
  benzene: 'Aromatic six-membered ring with alternating double bonds. Click a bond to fuse rings (naphthalene), or an atom to attach.',
  r3: 'Three-membered ring (cyclopropane / epoxide skeleton).', r4: 'Four-membered ring (cyclobutane, β-lactam skeleton).',
  r5: 'Five-membered ring (cyclopentane, furanose, pyrrolidine skeleton).', r6: 'Six-membered ring (cyclohexane, pyranose, piperidine skeleton).',
  r7: 'Seven-membered ring.', r8: 'Eight-membered ring.',
  'el-*': 'Polymer end (*): marks where a repeat unit continues. Each * has exactly one bond and MolCanvas draws a bracket with n across it.',
  ptable: 'Pick any element from the periodic table.',
  plus: 'Add a positive formal charge (e.g. ammonium N⁺, carbocation C⁺). The allowed number of bonds changes with charge.',
  minus: 'Add a negative formal charge (e.g. alkoxide O⁻, carbanion C⁻).',
  arrow: 'Reaction arrow: reactants → products. Drag to draw; it snaps to 15°. Shortcut A.',
  equilibrium: 'Equilibrium arrows for reversible reactions.', retro: 'Retrosynthesis arrow: target ⇒ precursors.',
  resonance: 'Double-headed arrow between resonance structures.',
  text: 'Add labels and conditions. Click an atom to type an element or a group such as OMe, CO2H, Boc or Ph. Shortcut T.',
};
for (const el of ['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I']) {
  const v: Record<string, string> = { C: '4 bonds', H: '1 bond', N: '3 bonds (4 as N⁺)', O: '2 bonds (3 as O⁺)', S: '2, 4 or 6 bonds', P: '3 or 5 bonds', F: '1 bond', Cl: '1 bond (up to 7 in oxyanions)', Br: '1 bond (up to 7 in oxyanions)', I: '1 bond (up to 7, e.g. IBX)' };
  TOOL_EXPLAIN[`el-${el}`] = `Valence: ${v[el]}. Click an atom to change it or empty space to add one. Missing bonds are filled with implicit hydrogens.`;
}

const toolButtons = new Map<string, HTMLButtonElement>();
const toolsEl = $('tools');
for (const def of TOOL_DEFS) {
  if (def === 'sep') { const s = document.createElement('div'); s.className = 'sep'; toolsEl.append(s); continue; }
  const b = document.createElement('button');
  b.className = `tool ${def.cls ?? ''}`;
  b.title = def.title;
  b.innerHTML = def.icon ? ICONS[def.icon] : def.label!;
  if (def.label && def.cls === 'el') b.style.color = atomColor(def.label, { theme: LIGHT, colorAtoms: true }) === LIGHT.fg ? '' : atomColor(def.label, { theme: LIGHT, colorAtoms: true });
  b.addEventListener('click', (e) => {
    if (def.key === 'ptable') { showPeriodicTable(e.currentTarget as HTMLElement); return; }
    editor.setTool(def.tool);
  });
  const name = def.title.replace(/\s*\(.*?\)/, '').split('.')[0];
  explainOnHover(b, name, TOOL_EXPLAIN[def.key] ?? def.title);
  toolButtons.set(def.key, b);
  toolsEl.append(b);
}

function toolKey(t: Tool): string {
  switch (t.kind) {
    case 'bond': return t.stereo !== 'none' ? t.stereo : `b${t.order}`;
    case 'ring': return t.aromatic ? 'benzene' : `r${t.size}`;
    case 'atom': return toolButtons.has(`el-${t.el}`) ? `el-${t.el}` : 'ptable';
    case 'charge': return t.delta > 0 ? 'plus' : 'minus';
    case 'arrow': return t.arrow === 'forward' ? 'arrow' : t.arrow;
    default: return t.kind;
  }
}
editor.onToolChange = (t) => {
  const k = toolKey(t);
  for (const [key, b] of toolButtons) {
    b.classList.toggle('active', key === k);
    if (key === 'ptable') b.style.color = '';
  }
  document.querySelectorAll('.card.active').forEach((c) => c.classList.remove('active'));
  if (t.kind === 'template') document.querySelector(`.card[data-idx="${PRESETS.indexOf(t.preset)}"]`)?.classList.add('active');
};

// ---------- periodic table ----------

const PT_ROWS = [
  'H . . . . . . . . . . . . . . . . He',
  'Li Be . . . . . . . . . . B C N O F Ne',
  'Na Mg . . . . . . . . . . Al Si P S Cl Ar',
  'K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr',
  'Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe',
  'Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn',
  'Fr Ra Ac Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og',
  '-',
  '. . . Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu .',
  '. . . Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr .',
];
function showPeriodicTable(anchor: HTMLElement) {
  const pop = $('ptable');
  if (!pop.hidden) { pop.hidden = true; return; }
  const grid = document.createElement('div');
  grid.className = 'ptable';
  for (const row of PT_ROWS) {
    if (row === '-') { const g = document.createElement('div'); g.className = 'gap'; grid.append(g); continue; }
    for (const sym of row.split(' ')) {
      if (sym === '.') { grid.append(document.createElement('span')); continue; }
      const b = document.createElement('button');
      b.textContent = sym;
      b.addEventListener('click', () => {
        pop.hidden = true;
        editor.setTool({ kind: 'atom', el: sym });
        toast(`${sym}: click an atom to change it, or click empty space.`);
      });
      grid.append(b);
    }
  }
  pop.replaceChildren(grid);
  pop.hidden = false;
  const r = anchor.getBoundingClientRect();
  pop.style.left = `${r.right + 8}px`;
  pop.style.top = `${Math.max(8, Math.min(window.innerHeight - pop.offsetHeight - 8, r.top - 120))}px`;
}
document.addEventListener('pointerdown', (e) => {
  const pop = $('ptable');
  if (!pop.hidden && !pop.contains(e.target as Node) && !(e.target as HTMLElement).closest('.tool')) pop.hidden = true;
  const menu = $('export-menu');
  if (!menu.hidden && !menu.contains(e.target as Node)) menu.hidden = true;
});

// ---------- library ----------

let activeCat: string | null = null;
const catsEl = $('lib-cats');
$('lib-search-icon').innerHTML = ICONS.search;
function renderCats() {
  catsEl.replaceChildren();
  for (const c of [null, ...CATEGORIES]) {
    const b = document.createElement('button');
    b.className = `chip${c === activeCat ? ' active' : ''}`;
    b.textContent = c ?? 'All';
    b.addEventListener('click', () => { activeCat = c; renderCats(); renderLibrary(); });
    catsEl.append(b);
  }
}

const thumbCache = new Map<string, string>();
function thumb(p: Preset): string {
  const key = `${p.smiles}|${isDark()}|${editor.colorAtoms}`;
  let s = thumbCache.get(key);
  if (!s) {
    try {
      let d = presetFrag(p);
      if (p.kind === 'group') {
        // show the attachment point as a generic R group rather than a polymer bracket
        d = JSON.parse(JSON.stringify(d)) as Doc;
        for (const a of d.atoms) if (a.el === '*') a.el = 'R';
      }
      s = standaloneSvg(d, { theme: isDark() ? DARK : LIGHT, colorAtoms: editor.colorAtoms, pad: 10, scale: 0.62 });
    } catch {
      s = '<span class="muted small">preview unavailable</span>';
    }
    thumbCache.set(key, s);
  }
  return s;
}

const observer = new IntersectionObserver((entries) => {
  for (const en of entries) {
    if (!en.isIntersecting) continue;
    const el = en.target as HTMLElement;
    const p = PRESETS[Number(el.dataset.idx)];
    el.querySelector('.thumb')!.innerHTML = thumb(p);
    observer.unobserve(el);
  }
}, { root: $('lib-grid'), rootMargin: '200px' });

function renderLibrary() {
  const grid = $('lib-grid');
  observer.disconnect();
  grid.replaceChildren();
  const q = ($('lib-filter') as HTMLInputElement).value.trim().toLowerCase();
  const matches = PRESETS.map((p, i) => ({ p, i })).filter(({ p }) =>
    (!activeCat || p.cat === activeCat) &&
    (!q || p.name.toLowerCase().includes(q) || (p.aka ?? '').toLowerCase().includes(q) || p.cat.toLowerCase().includes(q)));
  if (!matches.length) {
    const e = document.createElement('div');
    e.className = 'lib-empty';
    e.innerHTML = 'No presets match.<br/>Try “Add from text” on the right to load any SMILES or name.';
    grid.append(e);
    return;
  }
  let lastCat = '';
  for (const { p, i } of matches) {
    if (p.cat !== lastCat && !activeCat) {
      const t = document.createElement('div');
      t.className = 'lib-cat-title';
      t.textContent = p.cat;
      grid.append(t);
      lastCat = p.cat;
    }
    const card = document.createElement('div');
    card.className = 'card';
    card.dataset.idx = String(i);
    card.title = `${p.name}${p.aka ? ` · ${p.aka}` : ''}\n${p.smiles}\n\nClick to stamp, double-click to insert at centre.`;
    card.innerHTML = `<div class="thumb"></div><div class="name">${p.name}</div>`;
    const kindHelp = p.kind === 'group' ? 'Substituent: click an atom on the canvas to attach it (R marks the attachment point).'
      : p.kind === 'polymer' ? 'Polymer repeat unit: the brackets mark where the chain continues; n is the degree of polymerisation.'
      : 'Click to stamp it on the canvas, double-click to insert it.';
    explainOnHover(card, p.name, `${p.aka ? `${p.aka}. ` : ''}${kindHelp} SMILES ${p.smiles}`);
    card.addEventListener('click', () => {
      editor.setTool({ kind: 'template', preset: p, frag: presetFrag(p) });
      toast(p.kind === 'group' ? `${p.name}: click an atom to attach it.` : `${p.name}: click on the canvas to place it.`);
    });
    card.addEventListener('dblclick', () => {
      editor.insertFragment(JSON.parse(JSON.stringify(presetFrag(p))), false);
      editor.setTool({ kind: 'select' });
    });
    grid.append(card);
    observer.observe(card);
  }
}
$('lib-filter').addEventListener('input', renderLibrary);
renderCats();
renderLibrary();

// ---------- inspector ----------

const fmt = (n: number | undefined, d = 2) => (n === undefined || Number.isNaN(n) ? '–' : n.toFixed(d));
const formulaHtml = (f: string) => f.replace(/([A-Za-z])(\d+)/g, '$1<sub>$2</sub>').replace(/(\d*[+−])$/, '<sup>$1</sup>');
const COMP_COLORS: Record<string, string> = { C: '#8e8e93', H: '#c7c7cc', N: '#2f5bea', O: '#e0322b', S: '#c49a00', P: '#e07000', F: '#1e9e3a', Cl: '#34c759', Br: '#a52a2a', I: '#8a1fa0', Si: '#9a7752' };

let lastSmiles = '';
function updateInspector() {
  const useSel = editor.sel.size > 0 && editor.selectionDoc().atoms.length > 0;
  const doc = useSel ? editor.selectionDoc() : editor.doc;
  $('scope-title').textContent = useSel ? 'Properties · selection' : 'Properties';
  $('empty-hint').style.display = editor.doc.atoms.length || editor.doc.arrows.length || editor.doc.texts.length ? 'none' : '';
  undoBtn.disabled = !editor.canUndo();
  redoBtn.disabled = !editor.canRedo();
  let pr = null;
  try { pr = properties(doc); } catch { pr = null; }
  const props = $('props');
  if (!pr) {
    props.innerHTML = '<div class="muted small">Draw or select a structure to see its formula, mass and properties.</div>';
  } else {
    const comp = pr.composition.map((c) => `<span style="width:${c.pct}%;background:${COMP_COLORS[c.el] ?? '#aeaeb2'}"></span>`).join('');
    const legend = pr.composition.map((c) => `<span><i style="background:${COMP_COLORS[c.el] ?? '#aeaeb2'}"></i>${c.el} ${c.pct.toFixed(2)}%</span>`).join('');
    props.innerHTML = `
      <div class="formula">${formulaHtml(pr.formula)}${pr.isPolymer ? '<span class="muted small"> per repeat unit</span>' : ''}</div>
      <div class="kv">
        <span class="k">${pr.isPolymer ? 'Repeat unit mass (M₀)' : 'Molecular weight'}</span><span class="v">${fmt(pr.mw)} g/mol</span>
        <span class="k">Exact mass</span><span class="v">${fmt(pr.exactMass, 4)}</span>
        ${pr.isPolymer ? '' : `
        <span class="k">cLogP</span><span class="v">${fmt(pr.logP)}</span>
        <span class="k">TPSA</span><span class="v">${fmt(pr.tpsa, 1)} Å²</span>
        <span class="k">H-bond donors / acceptors</span><span class="v">${pr.hbd ?? '–'} / ${pr.hba ?? '–'}</span>
        <span class="k">Rotatable bonds</span><span class="v">${pr.rotatable ?? '–'}</span>
        <span class="k">Stereocentres</span><span class="v">${pr.stereocenters ?? '–'}</span>`}
        ${pr.fragments > 1 ? `<span class="k">Fragments</span><span class="v">${pr.fragments}</span>` : ''}
      </div>
      ${comp ? `<div class="comp">${comp}</div><div class="comp-legend">${legend}</div>` : ''}`;
  }
  updatePolymer(pr?.isPolymer ? pr.mw : null);
  let smiles = '';
  try { smiles = doc.atoms.length ? toSmiles(doc) : ''; } catch { smiles = ''; }
  $('smiles').textContent = smiles || '—';
  if (smiles !== lastSmiles) { $('name-out').textContent = ''; lastSmiles = smiles; }
}

function updatePolymer(m0: number | null) {
  const el = $('polymer');
  const stars = editor.doc.atoms.filter((a) => a.el === '*');
  const n = stars[0]?.sub ?? 'n';
  const nNum = parseFloat(n);
  el.innerHTML = `
    ${m0 !== null ? `
      <div class="field"><label for="dp">Degree of polymerisation n</label><input id="dp" value="${n}" /></div>
      <div class="kv" style="margin-top:8px">
        <span class="k">Mₙ ≈ n × M₀</span><span class="v">${Number.isFinite(nNum) ? `${(nNum * m0).toLocaleString(undefined, { maximumFractionDigits: 0 })} g/mol` : 'set n'}</span>
      </div>
      <div class="muted small" style="margin-top:4px">End groups are ignored in this estimate.</div>` :
      '<div class="muted small">Draw a monomer (or pick one from the library) and turn it into a repeat unit. Use the <b>*</b> tool to mark polymer ends.</div>'}
    <button class="mini wide" id="polymerize">Polymerise C=C monomer → repeat unit</button>`;
  const dp = el.querySelector<HTMLInputElement>('#dp');
  dp?.addEventListener('change', () => editor.setRepeatCount(dp.value.trim()));
  el.querySelector('#polymerize')!.addEventListener('click', () => {
    if (!editor.polymerize()) toast('No non-ring C=C double bond found to polymerise.');
  });
}

let inspectorTimer = 0;
let saveTimer = 0;
editor.onChange = () => {
  clearTimeout(inspectorTimer);
  inspectorTimer = window.setTimeout(updateInspector, 60);
  clearTimeout(saveTimer);
  saveTimer = window.setTimeout(() => {
    try { localStorage.setItem(STORAGE_KEY, JSON.stringify(editor.doc)); } catch { /* storage full or unavailable */ }
  }, 400);
};

$('copy-smiles').addEventListener('click', async () => {
  const s = $('smiles').textContent ?? '';
  if (!s || s === '—') return;
  try { await navigator.clipboard.writeText(s); toast('SMILES copied.'); } catch { toast('Copy failed: select the text and press ⌘C.'); }
});

// ---------- PubChem (optional, needs internet) ----------

const PUBCHEM = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';
async function nameToSmiles(name: string): Promise<string | null> {
  const r = await fetch(`${PUBCHEM}/compound/name/${encodeURIComponent(name)}/property/SMILES,IsomericSMILES/JSON`);
  if (!r.ok) return null;
  const j = await r.json();
  const p = j?.PropertyTable?.Properties?.[0];
  return p?.SMILES ?? p?.IsomericSMILES ?? null;
}
$('lookup-name').addEventListener('click', async () => {
  const s = $('smiles').textContent ?? '';
  const out = $('name-out');
  if (!s || s === '—') return;
  if (s.includes('*')) { out.textContent = 'Polymer repeat units are not in PubChem. Try the monomer instead.'; return; }
  out.textContent = 'Looking up…';
  try {
    const r = await fetch(`${PUBCHEM}/compound/smiles/property/IUPACName,Title/JSON`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
      body: `smiles=${encodeURIComponent(s)}`,
    });
    const j = await r.json();
    const p = j?.PropertyTable?.Properties?.[0];
    if (!p || (!p.Title && !p.IUPACName)) out.textContent = 'Not found in PubChem.';
    else out.innerHTML = `${p.Title ? `<b>${p.Title}</b><br/>` : ''}${p.IUPACName ?? ''}${p.CID ? `<br/>CID ${p.CID}` : ''}`;
  } catch {
    out.textContent = 'Could not reach PubChem (offline?).';
  }
});

$('import-btn').addEventListener('click', async () => {
  const ta = $('import-text') as HTMLTextAreaElement;
  const text = ta.value.trim();
  if (!text) return;
  try {
    editor.insertText(text);
    ta.value = '';
    return;
  } catch { /* not SMILES or MOL: try a name */ }
  toast(`Looking up “${text}” on PubChem…`);
  try {
    const smi = await nameToSmiles(text);
    if (!smi) { toast(`“${text}” is not valid SMILES and was not found on PubChem.`); return; }
    editor.insertText(smi);
    ta.value = '';
    toast(`Added ${text}.`);
  } catch {
    toast('Not valid SMILES/MOL, and PubChem could not be reached.');
  }
});

// ---------- files & export ----------

function newDrawing() {
  const hadContent = editor.doc.atoms.length > 0 || editor.doc.arrows.length > 0 || editor.doc.texts.length > 0;
  editor.clear();
  editor.fitView();
  if (hadContent) toast(`Canvas cleared. Press ${mod}Z to bring it back.`);
}

async function openFile() {
  const f = await openTextFile('.mcv,.json,.mol,.sdf,.smi,.smiles,.txt');
  if (!f) return;
  try {
    if (/\.(mcv|json)$/i.test(f.name)) {
      const d = JSON.parse(f.text) as Doc;
      if (!Array.isArray(d.atoms)) throw new Error('not a MolCanvas file');
      editor.load({ ...emptyDoc(), ...d });
    } else {
      editor.clear();
      editor.insertText(f.name.match(/\.(smi|smiles|txt)$/i) ? f.text.split('\n')[0] : f.text.split('$$$$')[0]);
      editor.fitView();
    }
    toast(`Opened ${f.name}`);
  } catch (e) {
    toast(`Could not open ${f.name}: ${(e as Error).message}`);
  }
}

async function saveDrawing() {
  const ok = await saveFile('structure.mcv', JSON.stringify(editor.doc, null, 1), { name: 'MolCanvas drawing', extensions: ['mcv'] });
  if (ok) toast('Saved.');
}

function exportDoc(): Doc {
  return editor.sel.size ? editor.selectionDoc() : editor.doc;
}

function exportSvg(background: string | null): string {
  return standaloneSvg(exportDoc(), { theme: LIGHT, colorAtoms: editor.colorAtoms, background, pad: 14 });
}

async function svgToPng(svgText: string, scale = 4): Promise<Blob> {
  const img = new Image();
  const url = URL.createObjectURL(new Blob([svgText], { type: 'image/svg+xml' }));
  await new Promise<void>((res, rej) => { img.onload = () => res(); img.onerror = () => rej(new Error('render failed')); img.src = url; });
  const canvas = document.createElement('canvas');
  canvas.width = Math.ceil(img.width * scale);
  canvas.height = Math.ceil(img.height * scale);
  const ctx = canvas.getContext('2d')!;
  ctx.scale(scale, scale);
  ctx.drawImage(img, 0, 0);
  URL.revokeObjectURL(url);
  return new Promise((res, rej) => canvas.toBlob((b) => (b ? res(b) : rej(new Error('PNG encode failed'))), 'image/png'));
}

async function copyImage() {
  if (!exportDoc().atoms.length && !exportDoc().arrows.length) { toast('Nothing to copy yet.'); return; }
  try {
    const svg = exportSvg('#ffffff');
    const png = svgToPng(svg);
    await navigator.clipboard.write([new ClipboardItem({ 'image/png': png })]);
    toast('Image copied: paste it into Word, PowerPoint or Keynote.');
  } catch {
    toast('Your system blocked image copying. Use Export → PNG instead.');
  }
}

function showExportMenu(anchor: HTMLElement) {
  const menu = $('export-menu');
  if (!menu.hidden) { menu.hidden = true; return; }
  const items: [string, string, () => Promise<void>][] = [
    ['PNG image', 'high-res, white background', async () => {
      const blob = await svgToPng(exportSvg('#ffffff'));
      if (await saveFile('structure.png', new Uint8Array(await blob.arrayBuffer()), { name: 'PNG image', extensions: ['png'] })) toast('PNG exported.');
    }],
    ['PNG (transparent)', 'for slides', async () => {
      const blob = await svgToPng(exportSvg(null));
      if (await saveFile('structure.png', new Uint8Array(await blob.arrayBuffer()), { name: 'PNG image', extensions: ['png'] })) toast('PNG exported.');
    }],
    ['SVG vector', 'for papers & Illustrator', async () => {
      if (await saveFile('structure.svg', exportSvg(null), { name: 'SVG image', extensions: ['svg'] })) toast('SVG exported.');
    }],
    ['MOL file', 'V2000, for other software', async () => {
      if (await saveFile('structure.mol', toMolfile(exportDoc()), { name: 'MDL Molfile', extensions: ['mol'] })) toast('MOL exported.');
    }],
    ['SMILES', '.smi text file', async () => {
      if (await saveFile('structure.smi', `${toSmiles(exportDoc())}\n`, { name: 'SMILES', extensions: ['smi'] })) toast('SMILES exported.');
    }],
  ];
  menu.replaceChildren(...items.map(([label, hint, fn]) => {
    const b = document.createElement('button');
    b.innerHTML = `${label}<span>${hint}</span>`;
    b.addEventListener('click', () => { menu.hidden = true; fn().catch((e: Error) => toast(`Export failed: ${e.message}`)); });
    return b;
  }));
  menu.hidden = false;
  const r = anchor.getBoundingClientRect();
  menu.style.left = `${r.left}px`;
  menu.style.top = `${r.bottom + 6}px`;
}

window.addEventListener('keydown', (e) => {
  const t = e.target as HTMLElement;
  if (t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA')) return;
  if (!(e.metaKey || e.ctrlKey)) return;
  const k = e.key.toLowerCase();
  if (k === 's') { e.preventDefault(); saveDrawing(); }
  else if (k === 'o') { e.preventDefault(); openFile(); }
  else if (k === 'n') { e.preventDefault(); newDrawing(); }
  else if (k === 'c' && e.shiftKey) { e.preventDefault(); copyImage(); }
});

// ---------- restore last session ----------

try {
  const saved = localStorage.getItem(STORAGE_KEY);
  if (saved) {
    const d = JSON.parse(saved) as Doc;
    if (Array.isArray(d.atoms)) editor.load({ ...emptyDoc(), ...d }, false);
  }
} catch { /* ignore corrupt or unavailable storage */ }

editor.setTool({ kind: 'bond', order: 1, stereo: 'none' });
updateInspector();
requestAnimationFrame(() => editor.fitView());
(window as unknown as { molcanvas: Editor }).molcanvas = editor; // handy for debugging from the console
