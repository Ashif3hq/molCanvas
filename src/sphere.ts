// Layer 3 of the solubility tools: experimental calibration of a polymer's Hansen solubility sphere.
//
// R₀ cannot be calculated from the structure: it comes from experiment. Given solvents the user has tested
// and scored as good (dissolves) or poor (does not), this fits the sphere centre (δD, δP, δH) and radius R₀
// with Hansen's SPHERE desirability function (C.M. Hansen, Hansen Solubility Parameters: A User's Handbook,
// 2nd ed., 2007, ch. 1):
//
//   Aᵢ = 1                 good solvent inside the sphere, or poor solvent outside
//   Aᵢ = exp(R₀ − Raᵢ)     good solvent outside
//   Aᵢ = exp(Raᵢ − R₀)     poor solvent inside
//   DATAFIT = (Π Aᵢ)^(1/n), maximised; among perfect fits the smallest sphere is taken.
//
// Once fitted, RED = Ra / R₀ for any other solvent: RED < 1 predicts dissolution, RED > 1 predicts none.
// The sphere belongs to the sample that was tested: molecular weight, crystallinity, temperature, cross-linking
// and the test itself (dissolved / swollen, concentration, time) all move it.
import { hansenDistance, HansenParameters } from './hansen';

export interface TestedSolvent {
  name: string;
  hsp: HansenParameters;
  good: boolean;
}

export interface SphereRow {
  name: string;
  good: boolean;
  ra: number;
  red: number;
  /** The fitted sphere disagrees with the experiment for this solvent. */
  misfit: boolean;
}

export interface Sphere {
  centre: HansenParameters;
  r0: number;
  /** Hansen's DATAFIT, 0–1; 1 means every solvent is on the right side of the sphere. */
  fit: number;
  rows: SphereRow[];
  nGood: number;
  nPoor: number;
  /** Things that make this sphere less trustworthy. */
  warnings: string[];
}

export type SphereResult = { ok: true; value: Sphere } | { ok: false; reason: string };

type Vec = [number, number, number, number]; // [2δD, δP, δH, R₀]: in this space Ra is a plain Euclidean distance

const toVec = (c: HansenParameters, r0: number): Vec => [2 * c.d, c.p, c.h, r0];
const centreOf = (x: Vec): HansenParameters => ({ d: x[0] / 2, p: x[1], h: x[2] });

function datafit(x: Vec, data: TestedSolvent[]): number {
  const c = centreOf(x), r0 = x[3];
  let logSum = 0;
  for (const s of data) {
    const ra = hansenDistance(c, s.hsp);
    if (s.good && ra > r0) logSum += r0 - ra;
    else if (!s.good && ra < r0) logSum += ra - r0;
  }
  return Math.exp(logSum / data.length);
}

/** Lower is better: −ln DATAFIT, a small pull towards the smallest sphere, and walls at δ < 0 and R₀ < 0.5. */
function objective(x: Vec, data: TestedSolvent[]): number {
  let wall = 0;
  for (let i = 0; i < 3; i++) if (x[i] < 0) wall += 1e3 * x[i] * x[i];
  if (x[3] < 0.5) wall += 1e3 * (0.5 - x[3]) ** 2;
  return -Math.log(datafit(x, data)) + 1e-4 * x[3] + wall;
}

function nelderMead(f: (x: Vec) => number, start: Vec, step: number[], iters = 4000): { x: Vec; fx: number } {
  const n = 4;
  let simplex: Vec[] = [start];
  for (let i = 0; i < n; i++) {
    const v = [...start] as Vec;
    v[i] += step[i];
    simplex.push(v);
  }
  let values = simplex.map(f);
  for (let it = 0; it < iters; it++) {
    const order = values.map((v, i) => [v, i] as const).sort((a, b) => a[0] - b[0]).map((p) => p[1]);
    simplex = order.map((i) => simplex[i]);
    values = order.map((i) => values[i]);
    if (Math.abs(values[n] - values[0]) < 1e-10 && it > 50) break;
    const centroid = [0, 0, 0, 0] as Vec;
    for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) centroid[k] += simplex[i][k] / n;
    const along = (t: number): Vec => centroid.map((c, k) => c + t * (simplex[n][k] - c)) as Vec;
    const xr = along(-1), fr = f(xr);
    if (fr < values[0]) {
      const xe = along(-2), fe = f(xe);
      if (fe < fr) { simplex[n] = xe; values[n] = fe; } else { simplex[n] = xr; values[n] = fr; }
    } else if (fr < values[n - 1]) {
      simplex[n] = xr; values[n] = fr;
    } else {
      const xc = fr < values[n] ? along(-0.5) : along(0.5), fc = f(xc);
      if (fc < Math.min(fr, values[n])) { simplex[n] = xc; values[n] = fc; }
      else {
        for (let i = 1; i <= n; i++) {
          simplex[i] = simplex[i].map((v, k) => simplex[0][k] + 0.5 * (v - simplex[0][k])) as Vec;
          values[i] = f(simplex[i]);
        }
      }
    }
  }
  let best = 0;
  for (let i = 1; i <= n; i++) if (values[i] < values[best]) best = i;
  return { x: simplex[best], fx: values[best] };
}

/**
 * Fit a Hansen solubility sphere to good/poor solvent results.
 * `guess` (for example the Layer 1 group-contribution estimate) is used as one extra starting point; it does not
 * bias the result.
 */
export function fitSphere(data: TestedSolvent[], guess?: HansenParameters): SphereResult {
  const good = data.filter((s) => s.good), poor = data.filter((s) => !s.good);
  if (!good.length) return { ok: false, reason: 'Mark at least one solvent as good (it dissolves the sample).' };
  if (!poor.length) return { ok: false, reason: 'Mark at least one solvent as poor: without non-solvents the sphere has no edge.' };

  const f = (x: Vec) => objective(x, data);
  const mean = (list: TestedSolvent[]): HansenParameters => ({
    d: list.reduce((a, s) => a + s.hsp.d, 0) / list.length,
    p: list.reduce((a, s) => a + s.hsp.p, 0) / list.length,
    h: list.reduce((a, s) => a + s.hsp.h, 0) / list.length,
  });
  const radiusFor = (c: HansenParameters) => Math.max(1, ...good.map((s) => hansenDistance(c, s.hsp))) + 0.1;
  const starts: HansenParameters[] = [mean(good), ...good.map((s) => s.hsp)];
  if (guess) starts.push(guess);

  let best: { x: Vec; fx: number } | null = null;
  for (const c of starts) {
    let r = nelderMead(f, toVec(c, radiusFor(c)), [2, 2, 2, 1]);
    r = nelderMead(f, r.x, [0.5, 0.5, 0.5, 0.3]); // restart to escape a collapsed simplex
    if (!best || r.fx < best.fx) best = r;
  }
  const x = best!.x;
  const centre = centreOf(x);
  // Report the radius to 0.01 and keep good solvents that sit on the surface inside it.
  const r0 = Math.max(0.5, x[3]);
  const rows: SphereRow[] = data.map((s) => {
    const ra = hansenDistance(centre, s.hsp);
    const inside = ra <= r0 + 1e-6;
    return { name: s.name, good: s.good, ra, red: ra / r0, misfit: inside !== s.good };
  });
  const fit = datafit(x, data);

  const warnings: string[] = [];
  if (good.length < 3) warnings.push(`only ${good.length} good solvent${good.length === 1 ? '' : 's'}: the centre is poorly defined`);
  if (poor.length < 3) warnings.push(`only ${poor.length} poor solvent${poor.length === 1 ? '' : 's'}: R₀ is poorly defined`);
  if (data.length < 10) warnings.push(`${data.length} solvents tested; Hansen recommends 20 or more chosen to cover δD, δP and δH`);
  const misfits = rows.filter((r) => r.misfit).length;
  if (misfits) warnings.push(`${misfits} solvent${misfits === 1 ? ' is' : 's are'} on the wrong side of the best sphere (check those results, or the sample may not follow a single sphere)`);
  // Is every direction bounded by a poor solvent? Check the poor solvents surround the centre in δP and δH.
  const around = (k: 'p' | 'h') => poor.some((s) => s.hsp[k] > centre[k]) && poor.some((s) => s.hsp[k] < centre[k]);
  for (const k of ['p', 'h'] as const) {
    if (!around(k)) warnings.push(`no poor solvent on one side in δ${k.toUpperCase()}: the sphere is open in that direction`);
  }

  return { ok: true, value: { centre, r0, fit, rows: rows.sort((a, b) => a.ra - b.ra), nGood: good.length, nPoor: poor.length, warnings } };
}
