// Layer 2 of the solubility tools: the Hansen comparison engine.
//
// This file knows nothing about drawn structures. It compares Hansen parameter sets that were obtained
// independently (from Layer 1, src/solubility.ts, from the solvent table below, from the user, or from a
// calibrated sphere in src/sphere.ts), and it works out the parameters of mixtures from their composition.
//
//   Ra² = 4(δD₁ − δD₂)² + (δP₁ − δP₂)² + (δH₁ − δH₂)²      (Hansen distance, MPa½)
//   RED = Ra / R₀                                           (only meaningful with an experimental R₀)
//
// Mixtures use the volume-fraction average δᵢ = Σ φⱼ δᵢ,ⱼ (Hansen, 2007, ch. 1), assuming ideal volume
// additivity. Mole or mass fractions are first converted to volume fractions with each component's molar
// volume (and molar mass, for mass fractions).

/** Hansen parameters in MPa½. */
export interface HansenParameters {
  d: number;
  p: number;
  h: number;
}

export function hansenDistance(polymer: HansenParameters, solvent: HansenParameters): number {
  return Math.sqrt(
    4 * (polymer.d - solvent.d) ** 2 +
    (polymer.p - solvent.p) ** 2 +
    (polymer.h - solvent.h) ** 2,
  );
}

/** Relative energy difference Ra / R₀. R₀ must come from experiment (see src/sphere.ts), not from structure. */
export function red(centre: HansenParameters, r0: number, solvent: HansenParameters): number {
  return hansenDistance(centre, solvent) / r0;
}

export const totalHansen = (x: HansenParameters) => Math.sqrt(x.d * x.d + x.p * x.p + x.h * x.h);

export interface Solvent {
  name: string;
  hsp: HansenParameters;
  /** Molar volume, cm³/mol (25 °C). */
  v?: number;
  /** Molar mass, g/mol. */
  m?: number;
  /** Where the values come from, e.g. "Hansen 2007", "drawn structure (group contribution)", "blend". */
  source: string;
}

const HANSEN_2007 = 'Hansen 2007';
const s = (name: string, d: number, p: number, h: number, v: number, m: number): Solvent =>
  ({ name, hsp: { d, p, h }, v, m, source: HANSEN_2007 });

/**
 * Experimental Hansen parameters of common solvents at 25 °C, MPa½, with molar volume (cm³/mol) and molar mass.
 * C.M. Hansen, Hansen Solubility Parameters: A User's Handbook, 2nd ed., CRC Press (2007), Table A.1.
 */
export const SOLVENTS: Solvent[] = [
  s('n-Hexane', 14.9, 0.0, 0.0, 131.6, 86.18),
  s('n-Heptane', 15.3, 0.0, 0.0, 147.4, 100.2),
  s('Cyclohexane', 16.8, 0.0, 0.2, 108.7, 84.16),
  s('Toluene', 18.0, 1.4, 2.0, 106.8, 92.14),
  s('Benzene', 18.4, 0.0, 2.0, 89.4, 78.11),
  s('p-Xylene', 17.8, 1.0, 3.1, 123.3, 106.2),
  s('Styrene', 18.6, 1.0, 4.1, 115.6, 104.2),
  s('Chlorobenzene', 19.0, 4.3, 2.0, 102.1, 112.6),
  s('Nitrobenzene', 20.0, 8.6, 4.1, 102.7, 123.1),
  s('Carbon tetrachloride', 17.8, 0.0, 0.6, 97.1, 153.8),
  s('Chloroform', 17.8, 3.1, 5.7, 80.7, 119.4),
  s('Dichloromethane', 18.2, 6.3, 6.1, 63.9, 84.93),
  s('1,2-Dichloroethane', 19.0, 7.4, 4.1, 79.4, 98.96),
  s('Trichloroethylene', 18.0, 3.1, 5.3, 90.2, 131.4),
  s('Diethyl ether', 14.5, 2.9, 5.1, 104.8, 74.12),
  s('Tetrahydrofuran', 16.8, 5.7, 8.0, 81.7, 72.11),
  s('1,4-Dioxane', 19.0, 1.8, 7.4, 85.7, 88.11),
  s('Acetone', 15.5, 10.4, 7.0, 74.0, 58.08),
  s('Methyl ethyl ketone', 16.0, 9.0, 5.1, 90.1, 72.11),
  s('Methyl isobutyl ketone', 15.3, 6.1, 4.1, 125.8, 100.2),
  s('Cyclohexanone', 17.8, 6.3, 5.1, 104.0, 98.14),
  s('Ethyl acetate', 15.8, 5.3, 7.2, 98.5, 88.11),
  s('n-Butyl acetate', 15.8, 3.7, 6.3, 132.5, 116.2),
  s('γ-Butyrolactone', 19.0, 16.6, 7.4, 76.8, 86.09),
  s('Propylene carbonate', 20.0, 18.0, 4.1, 85.0, 102.1),
  s('Acetonitrile', 15.3, 18.0, 6.1, 52.6, 41.05),
  s('Nitromethane', 15.8, 18.8, 5.1, 54.3, 61.04),
  s('Pyridine', 19.0, 8.8, 5.9, 80.9, 79.10),
  s('N,N-Dimethylformamide', 17.4, 13.7, 11.3, 77.0, 73.09),
  s('N,N-Dimethylacetamide', 16.8, 11.5, 10.2, 92.5, 87.12),
  s('N-Methyl-2-pyrrolidone', 18.0, 12.3, 7.2, 96.5, 99.13),
  s('Dimethyl sulfoxide', 18.4, 16.4, 10.2, 71.3, 78.13),
  s('m-Cresol', 18.0, 5.1, 12.9, 104.7, 108.1),
  s('Acetic acid', 14.5, 8.0, 13.5, 57.1, 60.05),
  s('1-Butanol', 16.0, 5.7, 15.8, 91.5, 74.12),
  s('2-Propanol', 15.8, 6.1, 16.4, 76.8, 60.10),
  s('1-Propanol', 16.0, 6.8, 17.4, 75.2, 60.10),
  s('Ethanol', 15.8, 8.8, 19.4, 58.5, 46.07),
  s('Methanol', 15.1, 12.3, 22.3, 40.7, 32.04),
  s('Ethylene glycol', 17.0, 11.0, 26.0, 55.8, 62.07),
  s('Glycerol', 17.4, 12.1, 29.3, 73.3, 92.09),
  s('Formamide', 17.2, 26.2, 19.0, 39.8, 45.04),
  s('Water', 15.5, 16.0, 42.3, 18.0, 18.02),
];

export interface Compared {
  solvent: Solvent;
  ra: number;
  /** Present only when an R₀ was given. */
  red?: number;
}

/** Ra (and RED, when R₀ is known) of each solvent against the target, closest first. */
export function compareSolvents(target: HansenParameters, solvents: Solvent[], r0?: number): Compared[] {
  return solvents
    .map((solvent) => {
      const ra = hansenDistance(target, solvent.hsp);
      return r0 && r0 > 0 ? { solvent, ra, red: ra / r0 } : { solvent, ra };
    })
    .sort((a, b) => a.ra - b.ra);
}

// ---------- mixtures ----------

export type CompositionBasis = 'volume' | 'mole' | 'mass';

export interface MixtureComponent {
  name: string;
  hsp: HansenParameters;
  /** Fraction or parts of this component on the chosen basis (any scale: they are normalised). */
  amount: number;
  /** Molar volume, cm³/mol: needed for mole and mass bases. */
  v?: number;
  /** Molar mass, g/mol: needed for the mass basis. */
  m?: number;
}

export interface Mixture {
  hsp: HansenParameters;
  /** Volume fraction of each component, in input order. */
  phi: number[];
  /** Molar volume of the mixture, cm³/mol, when every component has one (Σ xᵢVᵢ). */
  v?: number;
}

export type MixtureResult = { ok: true; value: Mixture } | { ok: false; reason: string };

/** Hansen parameters of a mixture from its composition, by volume-fraction averaging. */
export function mixtureHSP(components: MixtureComponent[], basis: CompositionBasis): MixtureResult {
  const no = (reason: string): MixtureResult => ({ ok: false, reason });
  if (components.length < 2) return no('A mixture needs at least two components.');
  for (const c of components) {
    if (!(c.amount >= 0) || !Number.isFinite(c.amount)) return no(`Enter a fraction of zero or more for ${c.name}.`);
    if (basis !== 'volume' && !(c.v && c.v > 0)) return no(`${c.name} has no molar volume, so ${basis} fractions cannot be converted to volume fractions.`);
    if (basis === 'mass' && !(c.m && c.m > 0)) return no(`${c.name} has no molar mass, so mass fractions cannot be converted.`);
  }
  // Volume of each component per unit of the basis.
  const vol = components.map((c) => basis === 'volume' ? c.amount : basis === 'mole' ? c.amount * c.v! : c.amount * c.v! / c.m!);
  const total = vol.reduce((a, b) => a + b, 0);
  if (!(total > 0)) return no('The fractions add up to zero.');
  const phi = vol.map((x) => x / total);
  const hsp = { d: 0, p: 0, h: 0 };
  components.forEach((c, i) => { hsp.d += phi[i] * c.hsp.d; hsp.p += phi[i] * c.hsp.p; hsp.h += phi[i] * c.hsp.h; });
  const value: Mixture = { hsp, phi };
  if (components.every((c) => c.v && c.v > 0)) {
    // Mole fractions from volume fractions: xᵢ ∝ φᵢ / Vᵢ; mixture molar volume Σ xᵢVᵢ = 1 / Σ(φᵢ/Vᵢ).
    value.v = 1 / components.reduce((s2, c, i) => s2 + phi[i] / c.v!, 0);
  }
  return { ok: true, value };
}
