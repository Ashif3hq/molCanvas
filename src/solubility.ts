// Solubility parameters by group contribution.
//
// Hansen (δD, δP, δH): Hoftyzer–Van Krevelen method,
//   δD = ΣFd / V,  δP = √(ΣFp²) / V,  δH = √(ΣEh / V)
//   group values from D.W. van Krevelen & K. te Nijenhuis, Properties of Polymers, 4th ed. (2009), Table 7.10.
// Hildebrand δ: Fedors method, δ = √(ΣEcoh / V)
//   group values from R.F. Fedors, Polym. Eng. Sci. 14 (1974) 147.
// V is the molar volume: the sum of Fedors volume increments, or M/ρ when the user gives a density.
// All δ values are in MPa½ (= (J/cm³)½).
//
// For a polymer repeat unit (exactly two * ends) the two * bonds are joined to each other, so a group
// split across the repeat-unit boundary (e.g. the ester in *OCCOC(=O)c1ccc(cc1)C(=O)*) is still seen as one group.
// Structures containing a group the tables do not cover are refused rather than estimated.
import { Molecule } from 'openchemlib';
import { Doc } from './model';
import { toOCL } from './chem';

type Hvk = [fd: number, fp: number, eh: number];
type Fed = [ecoh: number, v: number];

interface GroupDef {
  label: string;
  hvk: Hvk | Hvk[] | null; // an array means "sum of these HVK groups" (no single HVK entry exists)
  fedors: Fed | null;
}

// HVK vectors used for groups that HVK only lists in parts.
const HVK_CO: Hvk = [290, 770, 2000];
const HVK_COO: Hvk = [390, 490, 7000];
const HVK_O: Hvk = [100, 400, 3000];
const HVK_NH: Hvk = [160, 210, 3100];
const HVK_NH2: Hvk = [280, 0, 8400];
const HVK_N: Hvk = [20, 800, 5000];

const G: Record<string, GroupDef> = {
  CH3: { label: '–CH₃', hvk: [420, 0, 0], fedors: [4710, 33.5] },
  CH2: { label: '–CH₂–', hvk: [270, 0, 0], fedors: [4940, 16.1] },
  CH: { label: '>CH–', hvk: [80, 0, 0], fedors: [3430, -1.0] },
  C: { label: '>C<', hvk: [-70, 0, 0], fedors: [1470, -19.2] },
  dCH2: { label: '=CH₂', hvk: [400, 0, 0], fedors: [4310, 28.5] },
  dCH: { label: '=CH–', hvk: [200, 0, 0], fedors: [4310, 13.5] },
  dC: { label: '=C<', hvk: [70, 0, 0], fedors: [4310, -5.5] },
  tCH: { label: '≡CH', hvk: null, fedors: [3850, 27.4] },
  tC: { label: '–C≡', hvk: null, fedors: [7070, 6.5] },
  Ph1: { label: '–C₆H₅ (phenyl)', hvk: [1430, 110, 0], fedors: [31940, 71.4] },
  Ph2: { label: '–C₆H₄– (phenylene)', hvk: [1270, 110, 0], fedors: [31940, 52.4] },
  Ph3: { label: 'C₆H₃ (trisubst. benzene)', hvk: null, fedors: [31940, 33.4] },
  Ph4: { label: 'C₆H₂ (tetrasubst. benzene)', hvk: null, fedors: [31940, 14.4] },
  Ph5: { label: 'C₆H (pentasubst. benzene)', hvk: null, fedors: [31940, -4.6] },
  Ph6: { label: 'C₆ (hexasubst. benzene)', hvk: null, fedors: [31940, -23.6] },
  ring5: { label: 'ring (≥ 5 atoms)', hvk: [190, 0, 0], fedors: [1050, 16] },
  ring3: { label: 'ring (3–4 atoms)', hvk: [190, 0, 0], fedors: [3140, 18] },
  COOH: { label: '–COOH', hvk: [530, 420, 10000], fedors: [27630, 28.5] },
  COO: { label: '–COO– (ester)', hvk: HVK_COO, fedors: [18000, 18.0] },
  CO: { label: '>C=O (ketone)', hvk: HVK_CO, fedors: [17370, 10.8] },
  CHO: { label: '–CHO', hvk: [470, 800, 4500], fedors: [21350, 22.3] },
  O: { label: '–O– (ether)', hvk: HVK_O, fedors: [3350, 3.8] },
  OH: { label: '–OH', hvk: [210, 500, 20000], fedors: [29800, 10.0] },
  OH2: { label: '–OH (diol, gem or vicinal)', hvk: [210, 500, 20000], fedors: [21850, 13.0] },
  NH2: { label: '–NH₂', hvk: HVK_NH2, fedors: [12560, 19.2] },
  NH: { label: '>NH', hvk: HVK_NH, fedors: [8370, 4.5] },
  N: { label: '>N–', hvk: HVK_N, fedors: [4190, -9.0] },
  Nd: { label: '–N= (imine)', hvk: null, fedors: [11720, 5.0] },
  CN: { label: '–C≡N', hvk: [430, 1100, 2500], fedors: [25530, 24.0] },
  NO2: { label: '–NO₂ (aliphatic)', hvk: [500, 1070, 1500], fedors: [29310, 24.0] },
  NO2ar: { label: '–NO₂ (aromatic)', hvk: [500, 1070, 1500], fedors: [15360, 32.0] },
  S: { label: '–S–', hvk: [440, 0, 0], fedors: [14150, 12.0] },
  SH: { label: '–SH', hvk: null, fedors: [14440, 28.0] },
  F: { label: '–F', hvk: [220, 0, 0], fedors: [4190, 18.0] },
  F2: { label: '–F (2 on one C)', hvk: [220, 0, 0], fedors: [3560, 20.0] },
  F3: { label: '–F (≥ 3 on one C)', hvk: [220, 0, 0], fedors: [2300, 22.0] },
  Cl: { label: '–Cl', hvk: [450, 550, 400], fedors: [11550, 24.0] },
  Cl2: { label: '–Cl (2 on one C)', hvk: [450, 550, 400], fedors: [9630, 26.0] },
  Cl3: { label: '–Cl (≥ 3 on one C)', hvk: [450, 550, 400], fedors: [7530, 27.3] },
  Br: { label: '–Br', hvk: [550, 0, 0], fedors: [15490, 30.0] },
  Br2: { label: '–Br (2 on one C)', hvk: [550, 0, 0], fedors: [12350, 31.0] },
  Br3: { label: '–Br (≥ 3 on one C)', hvk: [550, 0, 0], fedors: [10670, 32.4] },
  I: { label: '–I', hvk: null, fedors: [19050, 31.5] },
  CONH2: { label: '–CONH₂', hvk: [HVK_CO, HVK_NH2], fedors: [41860, 17.5] },
  CONH: { label: '–CONH– (amide)', hvk: [HVK_CO, HVK_NH], fedors: [33490, 9.5] },
  CON: { label: '–CON< (amide)', hvk: [HVK_CO, HVK_N], fedors: [29510, -7.7] },
  OCOO: { label: '–OCOO– (carbonate)', hvk: [HVK_COO, HVK_O], fedors: [17580, 22.0] },
  NHCOO: { label: '–NHCOO– (urethane)', hvk: [HVK_NH, HVK_COO], fedors: [26370, 18.5] },
  Si: { label: '>Si<', hvk: null, fedors: [3390, 0] },
};

// Fedors: a halogen on a C=C or aromatic carbon has 20 % less Ecoh and 4 cm³/mol more volume.
const SP2_HALOGEN = { ecohFactor: 0.8, dv: 4.0 };

export interface SolubilityGroup {
  label: string;
  n: number;
}

export interface Solubility {
  groups: SolubilityGroup[];
  /** Molar volume used, cm³/mol, and where it came from. */
  v: number;
  vFedors: number;
  vSource: 'fedors' | 'density';
  /** Hildebrand δ by Fedors, MPa½. */
  hildebrand?: number;
  /** Hansen parameters by Hoftyzer–Van Krevelen, MPa½. */
  hansen?: { d: number; p: number; h: number; t: number };
  /** Why one method could not be applied while the other could. */
  hvkMissing?: string;
  fedorsMissing?: string;
  isPolymer: boolean;
}

export type SolubilityResult = { ok: true; value: Solubility } | { ok: false; reason: string };

interface Nb { j: number; order: number }

const SYM = (z: number) => Molecule.cAtomLabel[z];

/**
 * Group-contribution solubility parameters for a drawing (one molecule or one polymer repeat unit).
 * `density` (g/cm³) and `mass` (g/mol) are optional: given both, V = mass / density replaces the Fedors volume.
 */
export function solubility(doc: Doc, opts: { density?: number; mass?: number } = {}): SolubilityResult {
  const no = (reason: string): SolubilityResult => ({ ok: false, reason });
  if (!doc.atoms.length) return no('Nothing drawn.');
  const { mol } = toOCL(doc);
  mol.ensureHelperArrays(Molecule.cHelperRings);
  const n = mol.getAllAtoms();

  // Heavy-atom graph, with explicit H folded into hydrogen counts.
  const isH = (i: number) => mol.getAtomicNo(i) === 1 && !mol.getAtomCharge(i) && !mol.getAtomMass(i);
  const nh = new Array<number>(n).fill(0);
  const adj: Nb[][] = Array.from({ length: n }, () => []);
  for (let i = 0; i < n; i++) if (!isH(i)) nh[i] = mol.getImplicitHydrogens(i);
  for (let b = 0; b < mol.getAllBonds(); b++) {
    const a1 = mol.getBondAtom(0, b), a2 = mol.getBondAtom(1, b);
    if (isH(a1) && isH(a2)) return no('H₂ is not covered by the group tables.');
    if (isH(a1)) { nh[a2]++; continue; }
    if (isH(a2)) { nh[a1]++; continue; }
    const order = mol.getBondOrder(b);
    adj[a1].push({ j: a2, order });
    adj[a2].push({ j: a1, order });
  }
  const heavy: number[] = [];
  for (let i = 0; i < n; i++) if (!isH(i)) heavy.push(i);

  // One connected structure only.
  const comp = components(heavy, adj);
  if (comp > 1) return no('Draw or select a single molecule: mixtures and salts are not covered.');

  // Polymer ends: join the two * neighbours with a virtual single bond.
  const stars = heavy.filter((i) => mol.getAtomicNo(i) === 0);
  let isPolymer = false;
  const realEdges = heavy.reduce((s, i) => s + adj[i].length, 0) / 2;
  const realAtoms = heavy.length;
  if (stars.length) {
    if (stars.length !== 2 || stars.some((s) => adj[s].length !== 1)) {
      return no('Open * attachment points: solubility needs a whole molecule or a repeat unit with two * ends.');
    }
    const [s1, s2] = stars;
    const a1 = adj[s1][0].j, a2 = adj[s2][0].j;
    if (adj[s1][0].order !== 1 || adj[s2][0].order !== 1) return no('Polymer * ends must be single bonds.');
    for (const [s, a, other] of [[s1, a1, a2], [s2, a2, a1]]) {
      const k = adj[a].findIndex((e) => e.j === s);
      adj[a][k] = { j: other, order: 1 };
    }
    adj[s1] = [];
    adj[s2] = [];
    isPolymer = true;
  }
  const atoms = heavy.filter((i) => mol.getAtomicNo(i) !== 0);
  const z = (i: number) => mol.getAtomicNo(i);
  const el = (i: number) => SYM(z(i));
  const charge = (i: number) => mol.getAtomCharge(i);
  const degree = (i: number) => adj[i].length;
  const dbl = (i: number) => adj[i].filter((e) => e.order === 2);
  const isTerminalO = (i: number) => z(i) === 8 && degree(i) === 1 && nh[i] === 0;
  const aromatic = (i: number) => mol.isAromaticAtom(i);

  for (const i of atoms) {
    if (mol.getAtomRadical(i)) return no(`Radical on ${el(i)} is not covered.`);
    if (![6, 7, 8, 9, 14, 16, 17, 35, 53].includes(z(i))) return no(`${el(i)} is not in the group-contribution tables.`);
  }

  const counts = new Map<string, number>();
  const add = (key: string, k = 1) => counts.set(key, (counts.get(key) ?? 0) + k);
  const used = new Set<number>();
  const claim = (...ids: number[]) => {
    for (const id of ids) {
      if (used.has(id)) return false;
    }
    ids.forEach((id) => used.add(id));
    return true;
  };
  const sp2Halogen: string[] = [];

  // 1. Benzene rings (six-membered, all-carbon, aromatic). Fused aromatics and heteroaromatics are refused.
  const rings = mol.getRingSet();
  let benzenoid = 0, small = 0;
  const aromaticInBenzene = new Set<number>();
  for (let r = 0; r < rings.getSize(); r++) {
    const ra = rings.getRingAtoms(r);
    const size = rings.getRingSize(r);
    if (size <= 4) small++;
    if (size !== 6 || !ra.every((i) => z(i) === 6 && aromatic(i))) continue;
    if (ra.some((i) => aromaticInBenzene.has(i))) return no('Fused aromatic rings (e.g. naphthalene) are not covered.');
    const inRing = new Set(ra);
    const subs = ra.filter((i) => adj[i].some((e) => !inRing.has(e.j)) || nh[i] === 0).length;
    if (subs === 0) return no('Unsubstituted benzene is not in the tables (only phenyl and substituted benzene groups are).');
    ra.forEach((i) => aromaticInBenzene.add(i));
    if (!claim(...ra)) return no('Overlapping aromatic rings are not covered.');
    add(`Ph${subs}`);
    benzenoid++;
  }
  for (const i of atoms) {
    if (aromatic(i) && !aromaticInBenzene.has(i)) return no("Heteroaromatic and non-benzene aromatic rings are not covered.");
  }

  // 2. Carbonyl-based groups.
  for (const c of atoms) {
    if (z(c) !== 6 || used.has(c)) continue;
    const dO = adj[c].filter((e) => e.order === 2 && isTerminalO(e.j));
    if (!dO.length) continue;
    if (dO.length > 1 || degree(c) !== 3 - nh[c]) return no('Isocyanates, ketenes, CO₂ and other cumulated carbonyls are not covered.');
    const o = dO[0].j;
    const others = adj[c].filter((e) => e.j !== o);
    if (others.some((e) => e.order !== 1)) return no('Ketenes and cumulated C=O are not covered.');
    const nb = others.map((e) => e.j);
    const isC = (i: number) => z(i) === 6;
    const singleO = nb.filter((i) => z(i) === 8);
    const singleN = nb.filter((i) => z(i) === 7);
    if (nh[c] === 1) {
      if (nb.length === 1 && isC(nb[0])) { if (!claim(c, o)) return no('Overlapping groups.'); add('CHO'); continue; }
      return no('Formates, formamides and formaldehyde are not covered.');
    }
    if (nh[c] > 1) return no('Formaldehyde is not covered.');
    const ok = (fn: () => boolean, what: string) => { if (!fn()) throw new Error(what); };
    try {
      if (singleO.length === 2) {
        ok(() => singleO.every((x) => nh[x] === 0 && degree(x) === 2 && !charge(x)), 'Carbonic acid and its salts are not covered.');
        ok(() => claim(c, o, ...singleO), 'Anhydride-like groups are not covered.');
        add('OCOO');
      } else if (singleO.length === 1 && singleN.length === 1) {
        const x = singleO[0], nn = singleN[0];
        ok(() => nh[x] === 0 && degree(x) === 2 && nh[nn] === 1 && degree(nn) === 2 && !charge(nn), 'Only –NH–COO– urethanes are covered.');
        ok(() => claim(c, o, x, nn), 'Overlapping urethane groups are not covered.');
        add('NHCOO');
      } else if (singleO.length === 1 && isC(nb.find((i) => i !== singleO[0])!)) {
        const x = singleO[0];
        if (charge(x)) throw new Error('Carboxylate anions are not covered.');
        if (nh[x] === 1 && degree(x) === 1) { ok(() => claim(c, o, x), 'Overlapping groups.'); add('COOH'); }
        else if (nh[x] === 0 && degree(x) === 2 && adj[x].every((e) => z(e.j) === 6 && e.order === 1)) {
          ok(() => claim(c, o, x), 'Anhydrides are not covered.');
          add('COO');
        } else throw new Error('Peroxy acids and unusual esters are not covered.');
      } else if (singleN.length === 1 && singleO.length === 0 && isC(nb.find((i) => i !== singleN[0])!)) {
        const nn = singleN[0];
        if (charge(nn) || adj[nn].some((e) => e.order !== 1) || !adj[nn].every((e) => z(e.j) === 6)) throw new Error('Unusual amide nitrogen is not covered.');
        ok(() => claim(c, o, nn), 'Imides and ureas are not covered.');
        add(nh[nn] === 2 ? 'CONH2' : nh[nn] === 1 ? 'CONH' : 'CON');
      } else if (nb.length === 2 && nb.every(isC)) {
        ok(() => claim(c, o), 'Overlapping groups.');
        add('CO');
      } else {
        throw new Error('Acid halides, ureas and other carbonyl groups are not covered.');
      }
    } catch (e) {
      return no((e as Error).message);
    }
  }

  // 3. Nitriles and nitro groups.
  for (const i of atoms) {
    if (used.has(i)) continue;
    if (z(i) === 6) {
      const t = adj[i].find((e) => e.order === 3 && z(e.j) === 7 && degree(e.j) === 1);
      if (t) {
        if (degree(i) !== 2 || nh[i]) return no('HCN is not covered.');
        claim(i, t.j);
        add('CN');
      }
    }
    if (z(i) === 7 && charge(i) === 1 && degree(i) === 3) {
      const os = adj[i].filter((e) => isTerminalO(e.j));
      const cs = adj[i].filter((e) => z(e.j) === 6 && e.order === 1);
      if (os.length === 2 && cs.length === 1 && os.some((e) => e.order === 2) && os.some((e) => e.order === 1 && charge(e.j) === -1)) {
        claim(i, ...os.map((e) => e.j));
        add(aromatic(cs[0].j) ? 'NO2ar' : 'NO2');
      }
    }
  }

  const isHydroxyl = (o: number) => z(o) === 8 && nh[o] === 1 && degree(o) === 1 && !charge(o);
  const sp2C = (c: number) => aromatic(c) || dbl(c).some((e) => z(e.j) === 6);

  // 4. Remaining single atoms.
  for (const i of atoms) {
    if (used.has(i)) continue;
    if (charge(i)) return no(`Charged ${el(i)} is not covered (only neutral structures and nitro groups).`);
    const h = nh[i];
    const d = degree(i);
    const orders = adj[i].map((e) => e.order);
    const singles = orders.every((o) => o === 1);
    const nbEls = adj[i].map((e) => z(e.j));
    switch (z(i)) {
      case 6: {
        if (h + d !== 4 - orders.reduce((s, o) => s + o - 1, 0)) return no('Carbon with unusual valence.');
        const d2 = adj[i].filter((e) => e.order === 2), d3 = adj[i].filter((e) => e.order === 3);
        if (d2.some((e) => z(e.j) !== 6) || d3.some((e) => z(e.j) !== 6)) return no('C=N, C=S and similar double bonds are not covered.');
        if (d3.length === 1 && d2.length === 0) add(h === 1 ? 'tCH' : 'tC');
        else if (d2.length === 1 && d3.length === 0) add(h === 2 ? 'dCH2' : h === 1 ? 'dCH' : 'dC');
        else if (singles) {
          if (h === 4) return no('Methane is not covered.');
          add(['C', 'CH', 'CH2', 'CH3'][h]);
        } else return no('Allenes and cumulated double bonds are not covered.');
        break;
      }
      case 8: {
        if (!singles) return no('Unusual oxygen (C=O outside a carbonyl group) is not covered.');
        if (nbEls.some((x) => x === 8 || x === 7)) return no('Peroxides and N–O bonds are not covered.');
        if (h === 1 && d === 1) {
          const c = adj[i][0].j;
          if (z(c) !== 6) return no('Only C–OH hydroxyls are covered.');
          const near = [c, ...adj[c].map((e) => e.j).filter((x) => z(x) === 6)];
          const multi = near.some((x) => adj[x].some((e) => e.j !== i && isHydroxyl(e.j)));
          add(multi ? 'OH2' : 'OH');
        } else if (h === 0 && d === 2) add('O');
        else return no('Water and hydroxide are not covered.');
        break;
      }
      case 7: {
        if (nbEls.some((x) => x !== 6)) return no('N bonded to a heteroatom (hydrazines, N-oxides…) is not covered.');
        if (singles && h + d === 3 && h < 3) add(['N', 'NH', 'NH2'][h]);
        else if (orders.filter((o) => o === 2).length === 1 && h + d === 2 && h === 0) add('Nd');
        else return no('Unusual nitrogen (or ammonia) is not covered.');
        break;
      }
      case 16: {
        if (!singles || nbEls.some((x) => x !== 6)) return no('Only thioethers (–S–) and thiols (–SH) are covered, not sulfoxides, sulfones or disulfides.');
        if (h === 0 && d === 2) add('S');
        else if (h === 1 && d === 1) add('SH');
        else return no('Unusual sulfur is not covered.');
        break;
      }
      case 14: {
        if (h !== 0 || d !== 4 || !singles) return no('Only fully substituted >Si< is covered.');
        add('Si');
        break;
      }
      default: {
        // F, Cl, Br, I
        if (d !== 1 || z(adj[i][0].j) !== 6) return no(`${el(i)} must be bonded to carbon.`);
        const c = adj[i][0].j;
        const same = adj[c].filter((e) => z(e.j) === z(i)).length;
        const sym = el(i);
        let key = sym;
        if (same >= 2 && sym !== 'I') key = `${sym}${Math.min(3, same)}`;
        if (same >= 2 && sym === 'I') return no('Several iodines on one carbon are not covered.');
        add(key);
        if (sp2C(c)) sp2Halogen.push(key);
      }
    }
    used.add(i);
  }

  // 5. Rings: cyclomatic number of the real graph, minus benzene rings (already inside the phenyl values).
  const cyclomatic = realEdges - realAtoms + comp;
  const aliphaticRings = cyclomatic - benzenoid;
  const smallAliphatic = Math.min(small, aliphaticRings);
  if (aliphaticRings - smallAliphatic > 0) add('ring5', aliphaticRings - smallAliphatic);
  if (smallAliphatic > 0) add('ring3', smallAliphatic);

  // Sum both methods.
  let ecoh = 0, vF = 0, fedorsMissing: string | undefined;
  let fd = 0, fp2 = 0, eh = 0, hvkMissing: string | undefined;
  for (const [key, k] of counts) {
    const g = G[key];
    if (!g) return no('Group not covered.');
    if (g.fedors) { ecoh += k * g.fedors[0]; vF += k * g.fedors[1]; }
    else fedorsMissing ??= `${g.label} has no Fedors value`;
    if (g.hvk) {
      const parts = Array.isArray(g.hvk[0]) ? (g.hvk as Hvk[]) : [g.hvk as Hvk];
      for (const [a, b, c] of parts) {
        fd += k * a;
        eh += k * c;
        for (let r = 0; r < k; r++) fp2 += b * b;
      }
    } else hvkMissing ??= `${g.label} is not in the Hoftyzer–Van Krevelen table`;
  }
  for (const key of sp2Halogen) {
    const f = G[key].fedors;
    if (!f) continue;
    ecoh -= (1 - SP2_HALOGEN.ecohFactor) * f[0];
    vF += SP2_HALOGEN.dv;
  }
  if (fedorsMissing && hvkMissing) return no(`${fedorsMissing}, and ${hvkMissing}.`);
  if (fedorsMissing) return no(`${fedorsMissing}, so the molar volume cannot be estimated.`);
  if (vF <= 0) return no('The Fedors molar volume came out non-positive for this structure.');

  let v = vF, vSource: Solubility['vSource'] = 'fedors';
  if (opts.density && opts.density > 0 && opts.mass && opts.mass > 0) { v = opts.mass / opts.density; vSource = 'density'; }

  const groups = [...counts.entries()].map(([key, k]) => ({ label: G[key].label, n: k }));
  if (sp2Halogen.length) groups.push({ label: 'halogen on C=C / aryl (Fedors correction)', n: sp2Halogen.length });
  const value: Solubility = { groups, v, vFedors: vF, vSource, isPolymer, hildebrand: Math.sqrt(ecoh / v) };
  if (hvkMissing) value.hvkMissing = hvkMissing;
  else {
    const d = fd / v, p = Math.sqrt(fp2) / v, h = Math.sqrt(eh / v);
    value.hansen = { d, p, h, t: Math.sqrt(d * d + p * p + h * h) };
  }
  return { ok: true, value };
}

function components(atoms: number[], adj: Nb[][]): number {
  const seen = new Set<number>();
  let c = 0;
  for (const a of atoms) {
    if (seen.has(a)) continue;
    c++;
    const stack = [a];
    seen.add(a);
    while (stack.length) {
      const cur = stack.pop()!;
      for (const e of adj[cur]) if (!seen.has(e.j)) { seen.add(e.j); stack.push(e.j); }
    }
  }
  return c;
}
