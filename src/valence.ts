// Valence rules and plain-language explanations of atoms and bonds.
import { Atom, Bond, Doc, atomById, radicalElectrons, repeatUnitEnds, unpairedElectrons } from './model';

export const ELEMENT_NAMES: Record<string, string> = {
  H: 'Hydrogen', He: 'Helium', Li: 'Lithium', Be: 'Beryllium', B: 'Boron', C: 'Carbon', N: 'Nitrogen', O: 'Oxygen',
  F: 'Fluorine', Ne: 'Neon', Na: 'Sodium', Mg: 'Magnesium', Al: 'Aluminium', Si: 'Silicon', P: 'Phosphorus',
  S: 'Sulfur', Cl: 'Chlorine', Ar: 'Argon', K: 'Potassium', Ca: 'Calcium', Ti: 'Titanium', Cr: 'Chromium',
  Mn: 'Manganese', Fe: 'Iron', Co: 'Cobalt', Ni: 'Nickel', Cu: 'Copper', Zn: 'Zinc', Ga: 'Gallium', Ge: 'Germanium',
  As: 'Arsenic', Se: 'Selenium', Br: 'Bromine', Kr: 'Krypton', Rb: 'Rubidium', Sr: 'Strontium', Zr: 'Zirconium',
  Mo: 'Molybdenum', Ru: 'Ruthenium', Rh: 'Rhodium', Pd: 'Palladium', Ag: 'Silver', Cd: 'Cadmium', In: 'Indium',
  Sn: 'Tin', Sb: 'Antimony', Te: 'Tellurium', I: 'Iodine', Xe: 'Xenon', Cs: 'Caesium', Ba: 'Barium', W: 'Tungsten',
  Os: 'Osmium', Ir: 'Iridium', Pt: 'Platinum', Au: 'Gold', Hg: 'Mercury', Pb: 'Lead', Bi: 'Bismuth',
};

const VALENCE_ELECTRONS: Record<string, number> = { B: 3, C: 4, Si: 4, N: 5, P: 5, As: 5, O: 6, S: 6, Se: 6, F: 7, Cl: 7, Br: 7, I: 7 };

export const elementName = (el: string) => (el === '*' ? 'Polymer end' : ELEMENT_NAMES[el] ?? el);

const METALS = new Set(['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga',
  'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re',
  'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
  'Th', 'U', 'Fr', 'Ra', 'Ac']);

/**
 * Highest number of bonds (sum of bond orders) an atom may form, following the octet rule and
 * the usual expanded valences of period-3+ elements. Formal charge shifts the limit the way it does
 * for isoelectronic species (N⁺ behaves like C, O⁺ like N, C⁻ like N, B⁻ like C).
 */
export function maxValence(el: string, charge: number): number {
  switch (el) {
    case '*':
    case 'R': return 1;
    case 'H': return charge === 0 ? 1 : 0;
    case 'He': case 'Ne': case 'Ar': return 0;
    case 'B': return 3 - charge;
    case 'C': case 'Si': case 'Ge': return 4 - Math.abs(charge);
    case 'N': return 3 + charge;
    case 'O': return 2 + charge;
    case 'F': return 1 + charge;
    case 'P': case 'As': case 'Sb': return charge > 0 ? 4 : 5;
    case 'S': case 'Se': case 'Te': return 6;
    case 'Cl': case 'Br': case 'I': return 7;
    case 'Kr': case 'Xe': return 8;
  }
  if (METALS.has(el)) return 8; // coordination compounds: allow up to 8 bonds
  return 8;
}

/** Usual (lowest) valence, used to explain implicit hydrogens. */
export function typicalValence(el: string, charge: number): number {
  switch (el) {
    case 'P': case 'As': return 3 + charge;
    case 'S': case 'Se': return 2 + charge;
    case 'Cl': case 'Br': case 'I': return 1;
    default: return maxValence(el, charge);
  }
}

/**
 * Bond limit for a specific atom: the charge-adjusted valence minus the electrons held as radical dots.
 * A carbon radical (•CH3) forms 3 bonds, a carbene (:CH2) 2, an oxyl radical (RO•) 1, a chlorine atom (Cl•) none.
 */
export function valenceLimit(atom: Atom): number {
  return maxValence(atom.el, atom.charge) - radicalElectrons(atom.radical);
}

export function bondOrderSum(doc: Doc, atomId: number): number {
  let s = 0;
  for (const b of doc.bonds) if (b.a === atomId || b.b === atomId) s += b.order;
  return s;
}

export interface ValenceProblem {
  atom: Atom;
  used: number;
  max: number;
}

export function valenceProblems(doc: Doc): ValenceProblem[] {
  const out: ValenceProblem[] = [];
  for (const a of doc.atoms) {
    const used = bondOrderSum(doc, a.id);
    const max = valenceLimit(a);
    if (used > max) out.push({ atom: a, used, max });
  }
  return out;
}

const plural = (n: number, w: string) => `${n} ${w}${n === 1 ? '' : 's'}`;

/** A one-sentence reason for the user when a change is refused. */
export function explainProblem(p: ValenceProblem): string {
  const { atom, used, max } = p;
  if (atom.el === '*') return 'A polymer end (*) marks where the chain continues, so it can only have one bond.';
  const name = elementName(atom.el);
  const parts = [
    atom.charge ? `a ${atom.charge > 0 ? '+' : '−'}${Math.abs(atom.charge) > 1 ? Math.abs(atom.charge) : ''} charge` : '',
    atom.radical ? radicalPhrase(atom.radical) : '',
  ].filter(Boolean);
  const ch = parts.length ? ` with ${parts.join(' and ')}` : '';
  const why = atom.radical ? ' Each radical electron takes the place of one bond.' : '';
  if (max < 0) return `${name}${ch} is impossible: it has too few valence electrons for that.${why}`;
  if (max === 0) return `${name}${ch} cannot form bonds.${why}`;
  return `${name}${ch} can form at most ${plural(max, 'bond')} (valence ${max}). This change would give it ${used}.${why}`;
}

function radicalPhrase(r: Atom['radical']): string {
  return r === 'doublet' ? 'one unpaired electron (radical)'
    : r === 'triplet' ? 'two unpaired electrons (triplet diradical)'
    : 'a non-bonding electron pair (singlet carbene-type)';
}

/** Common name for a charged or radical atom, e.g. "tertiary carbocation", "radical anion", "nitrene". */
function speciesName(doc: Doc, atom: Atom, hCount: number): string {
  const bonds = doc.bonds.filter((b) => b.a === atom.id || b.b === atom.id);
  const used = bonds.reduce((s, b) => s + b.order, 0);
  const subst = used + hCount;
  const q = atom.charge, r = atom.radical;
  const carbonClass = () => {
    if (atom.el !== 'C' || subst !== 3 || bonds.some((b) => b.order !== 1)) return '';
    const nc = bonds.filter((b) => atomById(doc, b.a === atom.id ? b.b : b.a)!.el === 'C').length;
    return ['methyl', 'primary', 'secondary', 'tertiary'][nc] + ' ';
  };
  if (r === 'doublet') {
    if (q > 0) return 'radical cation';
    if (q < 0) return 'radical anion';
    if (atom.el === 'C') return `${carbonClass()}carbon radical`;
    if (atom.el === 'O') return subst === 1 ? 'oxyl radical (alkoxyl, hydroxyl or peroxyl)' : 'oxygen radical';
    if (atom.el === 'N') return subst === 2 ? 'aminyl radical' : 'nitrogen radical';
    if (subst === 0) return `${elementName(atom.el).toLowerCase()} atom (radical)`;
    return `${elementName(atom.el).toLowerCase()}-centred radical`;
  }
  if (r) {
    const spin = r === 'triplet' ? 'triplet' : 'singlet';
    if (!q && ['C', 'Si'].includes(atom.el) && subst === 2) return `${spin} ${atom.el === 'C' ? 'carbene' : 'silylene'}`;
    if (!q && atom.el === 'N' && subst === 1) return `${spin} nitrene`;
    return r === 'triplet' ? 'diradical' : '';
  }
  if (atom.el === 'C' && q > 0) return `${carbonClass()}carbocation`;
  if (atom.el === 'C' && q < 0) return `${carbonClass()}carbanion`;
  if (atom.el === 'N' && q === 1 && subst === 4) return bonds.some((b) => b.order > 1) ? 'iminium or nitro-type N⁺' : 'ammonium cation';
  if (atom.el === 'O' && q === 1 && subst === 3) return 'oxonium cation';
  if (atom.el === 'O' && q === -1 && subst === 1) return 'oxyanion (alkoxide, phenoxide, carboxylate O⁻…)';
  if (['F', 'Cl', 'Br', 'I'].includes(atom.el) && q === -1 && subst === 0) return `${{ F: 'fluoride', Cl: 'chloride', Br: 'bromide', I: 'iodide' }[atom.el]} anion`;
  if (q > 0) return 'cation';
  if (q < 0) return 'anion';
  return '';
}

/** Hybridisation from bond orders (VSEPR-style estimate for drawing purposes). */
function hybridization(doc: Doc, atom: Atom, hCount: number): string {
  const bonds = doc.bonds.filter((b) => b.a === atom.id || b.b === atom.id);
  const triples = bonds.filter((b) => b.order === 3).length;
  const doubles = bonds.filter((b) => b.order === 2).length;
  if (atom.el === 'H' || METALS.has(atom.el)) return '';
  if (triples || doubles >= 2) return 'sp';
  if (doubles === 1) return 'sp²';
  if (atom.el === 'C' && atom.radical === 'doublet' && bonds.length + hCount === 3) return 'sp² (near-planar radical)';
  if (atom.el === 'C' && atom.radical && bonds.length + hCount === 2) return atom.radical === 'triplet' ? 'sp (bent, triplet)' : 'sp² (singlet)';
  if (atom.el === 'C' && atom.charge !== 0 && bonds.length + hCount === 3) return atom.charge > 0 ? 'sp² (carbocation)' : 'sp³ (carbanion)';
  if (['C', 'N', 'O', 'Si'].includes(atom.el)) return 'sp³';
  return '';
}

export interface Explanation {
  title: string;
  detail: string;
  warn?: boolean;
}

export function explainAtom(doc: Doc, atom: Atom, hCount: number, aromaticRing: boolean): Explanation {
  const nb = doc.bonds.filter((b) => b.a === atom.id || b.b === atom.id).length;
  const used = bondOrderSum(doc, atom.id);
  const max = valenceLimit(atom);
  if (atom.el === '*') {
    return {
      title: 'Polymer end (*)',
      detail: nb <= 1 && repeatUnitEnds(doc).has(atom.id)
        ? 'Marks where the repeat unit continues into the next unit. MolCanvas draws a bracket across this bond and puts n on the last bracket. Click it with the Text tool to change n.'
        : nb <= 1
        ? 'An open bond or attachment point. A polymer repeat unit needs two * ends in the same molecule; add a second one and brackets with n appear.'
        : `A polymer end should have exactly one bond, but this one has ${nb}. Delete the extra bonds or change it to an element (hover and press C, N, O…).`,
      warn: nb > 1,
    };
  }
  const name = elementName(atom.el);
  const parts: string[] = [];
  const species = speciesName(doc, atom, hCount);
  if (species) parts.push(species);
  parts.push(`${plural(nb, 'bond')} (valence used ${used}/${Math.max(0, max)})`);
  if (hCount) parts.push(`${hCount} implicit H`);
  const lp = lonePairs(atom, used, hCount);
  if (lp > 0) parts.push(`${plural(lp, 'lone pair')}`);
  const unpaired = unpairedElectrons(atom.radical);
  if (unpaired) parts.push(`${plural(unpaired, 'unpaired electron')}${unpaired === 2 ? ' (triplet)' : ''}`);
  const hyb = hybridization(doc, atom, hCount);
  if (hyb) parts.push(aromaticRing && hyb === 'sp²' ? 'sp², aromatic ring' : hyb);
  if (atom.charge) parts.push(`formal charge ${atom.charge > 0 ? '+' : '−'}${Math.abs(atom.charge)}`);
  if (atom.iso) parts.push(`isotope ${atom.iso}${atom.el}`);
  const warn = used > max;
  const hint = warn
    ? `Too many bonds: ${explainProblem({ atom, used, max })}`
    : 'Press a letter (C, N, O, S…) to change the element, + or − for charge, . for a radical, Delete to remove.';
  return { title: `${name} (${atom.el})`, detail: `${parts.join(' · ')}. ${hint}`, warn };
}

/** Lone pairs from the octet count: valence electrons − charge − bonds − H − unpaired electrons, halved. */
export function lonePairs(atom: Atom, used: number, hCount: number): number {
  const ve = VALENCE_ELECTRONS[atom.el];
  if (ve === undefined) return 0;
  return Math.max(0, Math.floor((ve - atom.charge - used - hCount - unpairedElectrons(atom.radical)) / 2));
}

/** Recognise a few common functional groups so the info bar can name them. */
function bondGroupName(doc: Doc, b: Bond, A: Atom, B: Atom): string {
  const els = [A.el, B.el].sort().join('');
  if (b.order === 2 && els === 'CO') return 'carbonyl (C=O)';
  if (b.order === 2 && els === 'CC') return 'alkene (C=C)';
  if (b.order === 3 && els === 'CC') return 'alkyne (C≡C)';
  if (b.order === 3 && els === 'CN') return 'nitrile (C≡N)';
  if (b.order === 2 && els === 'CN') return 'imine (C=N)';
  if (b.order === 2 && els === 'NN') return 'azo (N=N)';
  if (b.order === 2 && els === 'OS') return 'sulfonyl/sulfinyl (S=O)';
  if (b.order === 2 && els === 'OP') return 'phosphoryl (P=O)';
  if (b.order === 1 && els === 'CO') return 'C–O single bond (alcohol, ether or ester)';
  if (b.order === 1 && els === 'CN') return 'C–N single bond (amine or amide)';
  if (els.includes('*')) return 'bond to a polymer end, crossed by the repeat-unit bracket';
  return '';
}

export function explainBond(doc: Doc, b: Bond, inRing: boolean): Explanation {
  const A = atomById(doc, b.a)!, B = atomById(doc, b.b)!;
  const kind = b.stereo === 'wedge' ? 'Wedge bond (towards you)'
    : b.stereo === 'hash' ? 'Hashed bond (away from you)'
    : b.stereo === 'wavy' && b.order === 2 ? 'Crossed double bond (E/Z unknown)'
    : b.stereo === 'wavy' ? 'Wavy bond (unknown stereochemistry)'
    : b.order === 3 ? 'Triple bond' : b.order === 2 ? 'Double bond' : 'Single bond';
  const sym = b.order === 3 ? '≡' : b.order === 2 ? '=' : '–';
  const group = bondGroupName(doc, b, A, B);
  const extras = [`${A.el}${sym}${B.el}`, group, inRing ? 'in a ring' : ''].filter(Boolean).join(' · ');
  return {
    title: kind,
    detail: `${extras}. Press 1, 2 or 3 to set the bond order, Delete to remove. With a bond tool, click to cycle single → double → triple.`,
  };
}
