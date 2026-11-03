// The solvent analyzer dialog: Layer 2 (compare a target with solvents) and Layer 3 (fit a solubility sphere
// from the user's good/poor results, then give RED). The target comes from Layer 1 (the drawn structure or a
// drawn mixture), from values the user types in (e.g. literature HSP and R₀), or from the fitted sphere.
import { compareSolvents, HansenParameters, mixtureHSP, Solvent, SOLVENTS } from './hansen';
import { fitSphere, Sphere } from './sphere';

export interface DrawnTarget {
  label: string;
  hsp: HansenParameters;
  /** Molar volume and mass, so the drawn structure can also be used as a solvent or blend component. */
  v?: number;
  m?: number;
  isPolymer: boolean;
}

type Mark = 'good' | 'poor';
type TargetMode = 'drawn' | 'sphere' | 'manual';

interface Saved {
  custom: Solvent[];
  marks: Record<string, Mark>;
  manual: { d: string; p: string; h: string; r0: string };
}

const STORE = 'molcanvas.solvents.v1';
const esc = (t: string) => t.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const f1 = (n: number) => (Number.isFinite(n) ? n.toFixed(1) : '–');
const f2 = (n: number) => (Number.isFinite(n) ? n.toFixed(2) : '–');
const num = (s: string) => { const v = parseFloat(s.replace(',', '.')); return Number.isFinite(v) ? v : NaN; };
const U = ' MPa<sup>½</sup>';

let saved: Saved = { custom: [], marks: {}, manual: { d: '', p: '', h: '', r0: '' } };
try {
  const raw = localStorage.getItem(STORE);
  if (raw) saved = { ...saved, ...JSON.parse(raw) };
} catch { /* storage unavailable or corrupt */ }
const persist = () => { try { localStorage.setItem(STORE, JSON.stringify(saved)); } catch { /* ignore */ } };

let drawn: DrawnTarget | null = null;
let sphere: Sphere | null = null;
let mode: TargetMode = 'drawn';
let dialog: HTMLDialogElement | null = null;
let notify: (msg: string) => void = () => {};

export function setNotifier(fn: (msg: string) => void) { notify = fn; }

/** Called by the inspector whenever Layer 1 (or the mixture calculation) has a new result for the drawing. */
export function setDrawnTarget(t: DrawnTarget | null) {
  drawn = t;
  if (dialog?.open) render();
}

const allSolvents = () => [...SOLVENTS, ...saved.custom];

function target(): { hsp: HansenParameters; r0?: number; label: string } | { error: string } {
  if (mode === 'sphere') {
    return sphere ? { hsp: sphere.centre, r0: sphere.r0, label: 'Fitted sphere (your experimental data)' } : { error: 'Fit a sphere first: mark good and poor solvents below.' };
  }
  if (mode === 'manual') {
    const m = saved.manual;
    const hsp = { d: num(m.d), p: num(m.p), h: num(m.h) };
    if (![hsp.d, hsp.p, hsp.h].every((x) => x >= 0)) return { error: 'Enter δD, δP and δH (MPa½) for the target.' };
    const r0 = num(m.r0);
    return { hsp, r0: r0 > 0 ? r0 : undefined, label: 'Your values' };
  }
  return drawn ? { hsp: drawn.hsp, label: `${drawn.label} (group-contribution estimate)` } : { error: 'Draw a structure with Hansen parameters, or enter values.' };
}

export function openAnalyzer() {
  if (!dialog) dialog = build();
  if (!drawn && mode === 'drawn') mode = sphere ? 'sphere' : 'manual';
  render();
  dialog.showModal();
}

function build(): HTMLDialogElement {
  const d = document.createElement('dialog');
  d.className = 'analyzer';
  d.setAttribute('aria-label', 'Solvent analyzer');
  document.body.appendChild(d);
  d.addEventListener('click', (e) => { if (e.target === d) d.close(); });
  return d;
}

function render() {
  const d = dialog!;
  const scrollTop = d.querySelector('.an-body')?.scrollTop ?? 0;
  const t = target();
  const rows = 'error' in t ? null : compareSolvents(t.hsp, allSolvents(), t.r0);
  const marked = Object.keys(saved.marks).filter((n) => allSolvents().some((s) => s.name === n));
  const m = saved.manual;
  const sel = (v: TargetMode, label: string, disabled = false) =>
    `<label class="${disabled ? 'muted' : ''}"><input type="radio" name="an-target" value="${v}" ${mode === v ? 'checked' : ''} ${disabled ? 'disabled' : ''}/> ${label}</label>`;

  d.innerHTML = `
    <header class="an-head">
      <h2>Solvent analyzer</h2>
      <button class="mini" data-act="close" aria-label="Close">Close</button>
    </header>
    <div class="an-body">
      <section>
        <h3>Target</h3>
        <div class="an-targets">
          ${sel('drawn', drawn ? `Drawn structure: ${esc(drawn.label)}` : 'Drawn structure (no Hansen estimate yet)', !drawn)}
          ${sel('sphere', sphere ? 'Fitted sphere from your solvent tests' : 'Fitted sphere (fit one below)', !sphere)}
          ${sel('manual', 'Enter values (e.g. literature δD, δP, δH and R₀)')}
        </div>
        ${mode === 'manual' ? `
        <div class="an-inputs">
          <label>δD <input data-manual="d" value="${esc(m.d)}" inputmode="decimal" /></label>
          <label>δP <input data-manual="p" value="${esc(m.p)}" inputmode="decimal" /></label>
          <label>δH <input data-manual="h" value="${esc(m.h)}" inputmode="decimal" /></label>
          <label>R₀ <input data-manual="r0" value="${esc(m.r0)}" placeholder="optional" inputmode="decimal" /></label>
        </div>` : ''}
        ${'error' in t ? `<p class="muted small">${esc(t.error)}</p>` : `
        <p class="small">δD ${f1(t.hsp.d)}, δP ${f1(t.hsp.p)}, δH ${f1(t.hsp.h)}${U}${t.r0 ? `, R₀ ${f1(t.r0)}${U}` : ''}</p>
        ${t.r0 ? '' : `<p class="muted small an-note">Only the Hansen distance Ra is shown. RED = Ra/R₀ needs the sphere radius R₀, which cannot be calculated from the structure: fit it from your solvent tests below, or enter a measured R₀.</p>`}`}
      </section>

      <section>
        <h3>Solvent compatibility <span class="muted">· Layer 2</span></h3>
        <p class="muted small">Ranked by Hansen distance Ra² = 4ΔδD² + ΔδP² + ΔδH². Mark the solvents you have tested to calibrate the sphere.</p>
        <div class="an-table-wrap">
          <table class="an-table">
            <thead><tr><th>Solvent</th><th>δD</th><th>δP</th><th>δH</th><th>Ra</th><th>RED</th><th>Your test</th><th></th></tr></thead>
            <tbody>
              ${(rows ?? allSolvents().map((solvent) => ({ solvent, ra: NaN, red: undefined as number | undefined }))).map((r) => {
                const mk = saved.marks[r.solvent.name];
                const cls = r.red === undefined ? '' : r.red < 1 ? 'in' : 'out';
                const custom = saved.custom.includes(r.solvent);
                return `<tr>
                  <td title="${esc(r.solvent.source)}">${esc(r.solvent.name)}${custom ? ` <span class="muted">· ${esc(r.solvent.source)}</span>` : ''}</td>
                  <td>${f1(r.solvent.hsp.d)}</td><td>${f1(r.solvent.hsp.p)}</td><td>${f1(r.solvent.hsp.h)}</td>
                  <td>${f1(r.ra)}</td>
                  <td class="${cls}">${r.red === undefined ? '–' : f2(r.red)}</td>
                  <td><select data-mark="${esc(r.solvent.name)}" aria-label="Test result for ${esc(r.solvent.name)}">
                    <option value="" ${!mk ? 'selected' : ''}>–</option>
                    <option value="good" ${mk === 'good' ? 'selected' : ''}>good</option>
                    <option value="poor" ${mk === 'poor' ? 'selected' : ''}>poor</option>
                  </select></td>
                  <td>${custom ? `<button class="an-x" data-remove="${esc(r.solvent.name)}" title="Remove">×</button>` : ''}</td>
                </tr>`;
              }).join('')}
            </tbody>
          </table>
        </div>
        <p class="muted small an-note">Solvent values: C.M. Hansen, <i>Hansen Solubility Parameters: A User's Handbook</i>, 2nd ed. (2007), 25 °C. RED &lt; 1 predicts the solvent dissolves the sample; RED &gt; 1 predicts it does not.</p>

        <details class="an-add">
          <summary>Add a solvent or a solvent blend</summary>
          <div class="an-inputs">
            <label>Name <input id="an-name" placeholder="e.g. Anisole" /></label>
            <label>δD <input id="an-d" inputmode="decimal" /></label>
            <label>δP <input id="an-p" inputmode="decimal" /></label>
            <label>δH <input id="an-h" inputmode="decimal" /></label>
            <label>V <input id="an-v" placeholder="cm³/mol" inputmode="decimal" /></label>
            <button class="mini" data-act="add">Add solvent</button>
          </div>
          <div class="an-inputs">
            <button class="mini" data-act="add-drawn" ${drawn && !drawn.isPolymer ? '' : 'disabled'}>Add the drawn molecule as a solvent</button>
            <span class="muted small">${drawn?.isPolymer ? 'A polymer repeat unit is not a solvent.' : 'Uses its group-contribution estimate.'}</span>
          </div>
          <div class="an-inputs">
            <label>Blend <select id="an-ba">${allSolvents().map((s) => `<option>${esc(s.name)}</option>`).join('')}</select></label>
            <label><input id="an-bf" value="50" inputmode="decimal" style="width:48px" /> % by</label>
            <select id="an-basis"><option value="volume">volume</option><option value="mass">mass</option><option value="mole">moles</option></select>
            <label>with <select id="an-bb">${allSolvents().map((s, i) => `<option ${i === 1 ? 'selected' : ''}>${esc(s.name)}</option>`).join('')}</select></label>
            <button class="mini" data-act="blend">Add blend</button>
          </div>
          <p class="muted small an-note">A blend's δ values are the volume-fraction average of its components (mass and mole fractions are converted with molar volume and molar mass). Two non-solvents can form a good blend.</p>
        </details>
      </section>

      <section>
        <h3>Polymer solubility prediction <span class="muted">· Layer 3</span></h3>
        <p class="muted small">Calibrate with experiment: dissolve your sample in a range of solvents, mark each as good or poor above (${marked.filter((n) => saved.marks[n] === 'good').length} good, ${marked.filter((n) => saved.marks[n] === 'poor').length} poor marked), then fit. The fit finds the sphere centre and R₀ that put good solvents inside and poor ones outside (Hansen's SPHERE method).</p>
        <div class="btn-row">
          <button class="mini primary" data-act="fit">Fit solubility sphere</button>
          <button class="mini" data-act="clear" ${marked.length ? '' : 'disabled'}>Clear test results</button>
        </div>
        ${sphere ? sphereHtml(sphere) : ''}
        <p class="muted small an-note">Group-contribution HSP are approximations (typically ±10 %). The real sphere also depends on molecular weight, crystallinity, temperature, cross-linking, morphology and the test itself (dissolved or only swollen, concentration, time), so it belongs to the sample you tested.</p>
      </section>
    </div>`;
  d.querySelector('.an-body')!.scrollTop = scrollTop;
  wire(d);
}

function sphereHtml(s: Sphere): string {
  const est = drawn ? Math.sqrt(4 * (drawn.hsp.d - s.centre.d) ** 2 + (drawn.hsp.p - s.centre.p) ** 2 + (drawn.hsp.h - s.centre.h) ** 2) : null;
  const misfits = s.rows.filter((r) => r.misfit);
  return `
    <div class="kv an-result">
      <span class="k">Centre δD, δP, δH</span><span class="v">${f1(s.centre.d)}, ${f1(s.centre.p)}, ${f1(s.centre.h)}${U}</span>
      <span class="k">Radius R₀</span><span class="v">${f1(s.r0)}${U}</span>
      <span class="k">Fit (DATAFIT, 1 = perfect)</span><span class="v">${s.fit.toFixed(3)}</span>
      <span class="k">Solvents used</span><span class="v">${s.nGood} good, ${s.nPoor} poor</span>
      ${est !== null ? `<span class="k">Ra from the structure estimate</span><span class="v">${f1(est)}${U}</span>` : ''}
    </div>
    ${misfits.length ? `<p class="small an-warn">Wrong side of the sphere: ${misfits.map((r) => `${esc(r.name)} (${r.good ? 'good' : 'poor'}, RED ${f2(r.red)})`).join(', ')}.</p>` : ''}
    ${s.warnings.length ? `<ul class="small an-warn">${s.warnings.map((w) => `<li>${esc(w)}</li>`).join('')}</ul>` : ''}
    ${mode !== 'sphere' ? '<button class="mini wide" data-act="use-sphere">Use the fitted sphere as the target (shows RED)</button>' : ''}`;
}

function wire(d: HTMLDialogElement) {
  const q = <T extends HTMLElement>(s: string) => d.querySelector<T>(s)!;
  for (const r of d.querySelectorAll<HTMLInputElement>('input[name="an-target"]')) {
    r.addEventListener('change', () => { mode = r.value as TargetMode; render(); });
  }
  for (const inp of d.querySelectorAll<HTMLInputElement>('input[data-manual]')) {
    inp.addEventListener('change', () => { saved.manual[inp.dataset.manual as keyof Saved['manual']] = inp.value.trim(); persist(); render(); });
  }
  for (const s of d.querySelectorAll<HTMLSelectElement>('select[data-mark]')) {
    s.addEventListener('change', () => {
      const name = s.dataset.mark!;
      if (s.value) saved.marks[name] = s.value as Mark; else delete saved.marks[name];
      persist();
      render();
    });
  }
  for (const b of d.querySelectorAll<HTMLButtonElement>('button[data-remove]')) {
    b.addEventListener('click', () => {
      const name = b.dataset.remove!;
      saved.custom = saved.custom.filter((s) => s.name !== name);
      delete saved.marks[name];
      persist();
      render();
    });
  }
  const act = (name: string, fn: () => void) => d.querySelector(`[data-act="${name}"]`)?.addEventListener('click', fn);
  act('close', () => d.close());
  act('add', () => {
    const name = q<HTMLInputElement>('#an-name').value.trim();
    const hsp = { d: num(q<HTMLInputElement>('#an-d').value), p: num(q<HTMLInputElement>('#an-p').value), h: num(q<HTMLInputElement>('#an-h').value) };
    const v = num(q<HTMLInputElement>('#an-v').value);
    if (!name) return notify('Give the solvent a name.');
    if (![hsp.d, hsp.p, hsp.h].every((x) => x >= 0)) return notify('Enter δD, δP and δH in MPa½ (zero or more).');
    addCustom({ name, hsp, v: v > 0 ? v : undefined, source: 'your values' });
  });
  act('add-drawn', () => {
    if (!drawn) return;
    addCustom({ name: drawn.label, hsp: drawn.hsp, v: drawn.v, m: drawn.m, source: 'drawn, group contribution' });
  });
  act('blend', () => {
    const a = allSolvents().find((s) => s.name === q<HTMLSelectElement>('#an-ba').value)!;
    const b = allSolvents().find((s) => s.name === q<HTMLSelectElement>('#an-bb').value)!;
    const pct = num(q<HTMLInputElement>('#an-bf').value);
    const basis = q<HTMLSelectElement>('#an-basis').value as 'volume' | 'mass' | 'mole';
    if (a === b) return notify('Pick two different solvents for a blend.');
    if (!(pct > 0 && pct < 100)) return notify('The first solvent’s share must be between 0 and 100 %.');
    const r = mixtureHSP([{ name: a.name, hsp: a.hsp, v: a.v, m: a.m, amount: pct }, { name: b.name, hsp: b.hsp, v: b.v, m: b.m, amount: 100 - pct }], basis);
    if (!r.ok) return notify(r.reason);
    const unit = basis === 'volume' ? 'v/v' : basis === 'mass' ? 'w/w' : 'mol/mol';
    addCustom({ name: `${a.name}/${b.name} ${+pct.toFixed(1)}:${+(100 - pct).toFixed(1)} ${unit}`, hsp: r.value.hsp, v: r.value.v, source: `blend, ${(100 * r.value.phi[0]).toFixed(0)} vol % ${a.name}` });
  });
  act('fit', () => {
    const data = allSolvents().filter((s) => saved.marks[s.name]).map((s) => ({ name: s.name, hsp: s.hsp, good: saved.marks[s.name] === 'good' }));
    const r = fitSphere(data, drawn?.hsp);
    if (!r.ok) return notify(r.reason);
    sphere = r.value;
    mode = 'sphere';
    render();
  });
  act('clear', () => { saved.marks = {}; sphere = null; if (mode === 'sphere') mode = drawn ? 'drawn' : 'manual'; persist(); render(); });
  act('use-sphere', () => { mode = 'sphere'; render(); });
}

function addCustom(s: Solvent) {
  if (allSolvents().some((x) => x.name === s.name)) return notify(`“${s.name}” is already in the list.`);
  saved.custom.push(s);
  persist();
  render();
}
