// Checks ChemDraw CDXML import and export: round trips of structures, stereo, polymers and reaction
// schemes, plus a file written the way ChemDraw 21 writes them (abbreviations, aromatic bonds,
// superseded arrows, SRU brackets). Pass an output directory to also write every exported file there,
// e.g. for validation with RDKit (Chem.MolsFromCDXMLFile) or Open Babel.
import { mkdirSync, writeFileSync } from 'node:fs';
import { fromCdxml, parseXml, toCdxml } from '../src/cdxml';
import { Molecule } from 'openchemlib';
import { fromSmiles, toOCL, toSmiles } from '../src/chem';
import { ArrowKind, BOND, Doc, fragments, mergeInto, subDoc } from '../src/model';

const outDir = process.argv[2];
if (outDir) mkdirSync(outDir, { recursive: true });

const CASES: Record<string, string> = {
  ethanol: 'CCO',
  alanine: 'C[C@H](N)C(=O)O',
  ibuprofen: 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O',
  nitrobenzene: 'c1ccccc1[N+](=O)[O-]',
  tetramethylammonium: 'C[N+](C)(C)C',
  transStilbene: 'C(/c1ccccc1)=C\\c1ccccc1',
  cisButene: 'C/C=C\\C',
  deuteriochloroform: '[2H]C(Cl)(Cl)Cl',
  carbon13: '[13CH3]O',
  pyridine: 'c1ccncc1',
  pyrrole: 'c1cc[nH]c1',
  sulfone: 'CS(=O)(=O)c1ccccc1',
  phosphate: 'COP(=O)(OC)OC',
  acetonitrile: 'CC#N',
  pmma: '*CC(*)(C)C(=O)OC',
  pet: '*OCCOC(=O)c1ccc(cc1)C(=O)*',
  nylon66: '*NCCCCCCNC(=O)CCCCC(=O)*',
};

let bad = 0;
const fail = (name: string, msg: string) => { bad++; console.error('FAIL', name, msg); };
const idcodes = (d: Doc) => fragments(d).map((c) => {
  const { mol } = toOCL(subDoc(d, new Set(c)));
  mol.ensureHelperArrays(Molecule.cHelperParities);
  return mol.getIDCode();
}).sort().join(' ');
const smilesSet = (d: Doc) => fragments(d).map((c) => toSmiles(subDoc(d, new Set(c)))).sort().join('.');

for (const [name, smi] of Object.entries(CASES)) {
  const doc = fromSmiles(smi);
  const xml = toCdxml(doc);
  try { parseXml(xml); } catch (e) { fail(name, `exported XML does not parse: ${(e as Error).message}`); continue; }
  const { doc: back, notes } = fromCdxml(xml);
  if (toSmiles(back) !== toSmiles(doc)) fail(name, `round trip gave ${toSmiles(back)}, expected ${toSmiles(doc)}`);
  if (notes.length) fail(name, `unexpected import notes: ${notes.join(' ')}`);
  if (smi.includes('*') && !xml.includes('BracketUsage="SRU"')) fail(name, 'repeat unit was not written as an SRU bracket');
  const lens = back.bonds.map((b) => {
    const A = back.atoms.find((a) => a.id === b.a)!, B = back.atoms.find((a) => a.id === b.b)!;
    return Math.hypot(A.x - B.x, A.y - B.y);
  });
  if (lens.some((l) => Math.abs(l - BOND) > 1)) fail(name, 'bond lengths not restored to the standard length');
  console.log('ok  ', name.padEnd(20), toSmiles(back));
  if (outDir) writeFileSync(`${outDir}/${name}.cdxml`, xml);
}

// Unknown E/Z (crossed double bond) survives the trip and does not gain a configuration.
{
  const doc = fromSmiles('CC=CC');
  doc.bonds.find((b) => b.order === 2)!.stereo = 'wavy';
  const xml = toCdxml(doc);
  const back = fromCdxml(xml).doc;
  if (!xml.includes('Display="Wavy"') || back.bonds.find((b) => b.order === 2)?.stereo !== 'wavy' || /[/\\]/.test(toSmiles(back))) fail('either-butene', 'unknown double-bond geometry was lost');
  else console.log('ok  ', 'either-butene'.padEnd(20), toSmiles(back));
  if (outDir) writeFileSync(`${outDir}/either-butene.cdxml`, xml);
}

// A reaction scheme: reactants, arrow, conditions above the arrow, product; and every arrow kind.
{
  const doc = fromSmiles('CC(=O)O');
  mergeInto(doc, fromSmiles('CCO'), BOND * 3, 0);
  mergeInto(doc, fromSmiles('CC(=O)OCC'), BOND * 12, 0);
  const kinds: ArrowKind[] = ['forward', 'equilibrium', 'retro', 'resonance'];
  doc.arrows.push({ id: doc.nextId++, x1: BOND * 5, y1: 0, x2: BOND * 9, y2: 0, kind: 'forward' });
  kinds.slice(1).forEach((kind, i) => doc.arrows.push({ id: doc.nextId++, x1: 0, y1: BOND * (5 + 2 * i), x2: BOND * 3, y2: BOND * (5 + 2 * i), kind }));
  doc.texts.push({ id: doc.nextId++, x: BOND * 5.5, y: -BOND * 0.4, text: 'H2SO4' });
  doc.texts.push({ id: doc.nextId++, x: BOND * 1.4, y: 0, text: '+' });
  const xml = toCdxml(doc);
  const { doc: back, notes } = fromCdxml(xml);
  const step = parseXml(xml).children[0].children.find((c) => c.tag === 'page')!.children.find((c) => c.tag === 'scheme')?.children[0];
  if (!step || step.attrs.ReactionStepReactants.split(' ').length !== 2 || step.attrs.ReactionStepProducts.split(' ').length !== 1 || !step.attrs.ReactionStepObjectsAboveArrow) {
    fail('reaction', `scheme step wrong: ${JSON.stringify(step?.attrs)}`);
  }
  if (back.arrows.map((a) => a.kind).join() !== doc.arrows.map((a) => a.kind).join()) fail('reaction', 'arrow kinds changed');
  const dir = (d: Doc) => d.arrows.map((a) => Math.sign(a.x2 - a.x1)).join();
  if (dir(back) !== dir(doc)) fail('reaction', 'arrow directions changed');
  if (back.texts.map((t) => t.text).join() !== 'H2SO4,+') fail('reaction', `texts: ${back.texts.map((t) => t.text)}`);
  if (smilesSet(back) !== smilesSet(doc)) fail('reaction', `structures: ${smilesSet(back)}`);
  if (notes.length) fail('reaction', notes.join(' '));
  console.log('ok  ', 'reaction'.padEnd(20), smilesSet(back), back.arrows.map((a) => a.kind).join(' '));
  if (outDir) writeFileSync(`${outDir}/reaction.cdxml`, xml);
}

// Written the way ChemDraw 21 writes files: an OMe abbreviation (nested fragment with an external
// connection point), aromatic 1.5 bonds, an NH2 label, a hashed wedge stored as WedgedHashEnd, an old-style arrow
// superseded by an <arrow>, a reaction +, a free-standing "THF" label and a polyethylene SRU.
const CHEMDRAW = `<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
<CDXML CreationProgram="ChemDraw 21.0.0.28" BondLength="30" LabelFont="21" LabelSize="10"><colortable><color r="1" g="1" b="1"/><color r="0" g="0" b="0"/></colortable><fonttable><font id="21" charset="x-mac-roman" name="Helvetica"/></fonttable><page id="1" BoundingBox="0 0 540 720">
<fragment id="2">
<n id="3" p="100 100" AS="N"/><n id="4" p="125.98 85" AS="N"/><n id="5" p="151.96 100" AS="N"/><n id="6" p="151.96 130" AS="N"/><n id="7" p="125.98 145" AS="N"/><n id="8" p="100 130" AS="N"/>
<n id="9" p="177.94 85" NodeType="Nickname" NeedsClean="yes" AS="N"><fragment id="10"><n id="11" p="177.94 85" Element="8" NumHydrogens="0"/><n id="12" p="203.92 100"/><n id="13" p="151.96 100" NodeType="ExternalConnectionPoint" ExternalConnectionNum="1"/><b id="14" B="11" E="12"/><b id="15" B="13" E="11"/></fragment><t p="174.05 88.68"><s font="21" size="10" color="0" face="96">OMe</s></t></n>
<n id="16" p="74.02 85" AS="N"/><n id="17" p="48.04 100" Element="7" NumHydrogens="2" AS="N"><t p="44.43 103.68"><s font="21" size="10" color="0" face="96">NH2</s></t></n><n id="18" p="74.02 55" AS="N"/>
<b id="20" B="3" E="4" Order="1.5"/><b id="21" B="4" E="5" Order="1.5"/><b id="22" B="5" E="6" Order="1.5"/><b id="23" B="6" E="7" Order="1.5"/><b id="24" B="7" E="8" Order="1.5"/><b id="25" B="8" E="3" Order="1.5"/>
<b id="26" B="5" E="9"/><b id="27" B="17" E="16" Display="WedgedHashEnd"/><b id="28" B="16" E="3"/><b id="29" B="16" E="18"/>
</fragment>
<fragment id="40"><n id="41" p="300 100"/><n id="42" p="325.98 85"/><n id="43" p="351.96 100"/><n id="44" p="377.94 85"/>
<b id="45" B="41" E="42"/><b id="46" B="42" E="43"/><b id="47" B="43" E="44"/></fragment>
<graphic id="50" BoundingBox="300 90 300 110" GraphicType="Bracket" BracketType="Square"/>
<graphic id="51" BoundingBox="365 95 365 75" GraphicType="Bracket" BracketType="Square" BracketUsage="SRU"/>
<bracketedgroup id="52" BracketedObjectIDs="42 43" BracketUsage="SRU" SRULabel="n"><bracketattachment id="53" GraphicID="50"><crossingbond id="54" BondID="45" InnerAtomID="42"/></bracketattachment><bracketattachment id="55" GraphicID="51"><crossingbond id="56" BondID="47" InnerAtomID="43"/></bracketattachment></bracketedgroup>
<graphic id="60" SupersededBy="61" BoundingBox="280 200 220 200" GraphicType="Line" ArrowType="FullHead" HeadSize="1000"/>
<arrow id="61" BoundingBox="220 196 280 204" FillType="None" ArrowheadHead="Full" ArrowheadType="Solid" HeadSize="1000" Head3D="280 200 0" Tail3D="220 200 0"/>
<graphic id="62" BoundingBox="200 200 200 207.5" GraphicType="Symbol" SymbolType="Plus"/>
<fragment id="70"><n id="71" p="240 190" NodeType="Unspecified"><t p="236 194"><s font="21" size="10" color="0">THF</s></t></n></fragment>
</page></CDXML>`;
{
  const { doc, notes } = fromCdxml(CHEMDRAW);
  // compared by OpenChemLib ID code, which is canonical and includes stereo (checked against RDKit's reading)
  const got = idcodes(doc);
  const want = idcodes(fromSmiles('COc1cccc([C@@H](C)N)c1.*CC*'));
  if (got !== want) fail('chemdraw-file', `structures ${smilesSet(doc)}, expected COc1cccc([C@@H](C)N)c1 and *CC*`);
  if (doc.arrows.length !== 1 || doc.arrows[0].kind !== 'forward' || doc.arrows[0].x2 <= doc.arrows[0].x1) fail('chemdraw-file', 'arrow (head on the right) not read once');
  if (doc.texts.map((t) => t.text).sort().join() !== '+,THF') fail('chemdraw-file', `texts: ${doc.texts.map((t) => t.text)}`);
  if (!doc.atoms.some((a) => a.el === '*' && a.sub === 'n')) fail('chemdraw-file', 'SRU did not become a repeat unit');
  if (notes.length) fail('chemdraw-file', `unexpected notes: ${notes.join(' ')}`);
  console.log('ok  ', 'chemdraw-file'.padEnd(20), smilesSet(doc));
}

// Things MolCanvas cannot represent are reported, never silently changed.
{
  const xml = `<CDXML><page><fragment><n id="1" p="0 0"/><n id="2" p="14.4 0" NodeType="GenericNickname" GenericNickname="R"/><b id="3" B="1" E="2"/><n id="4" p="28.8 0"/><b id="5" B="1" E="4" Order="hydrogen"/></fragment>
<arrow id="6" Head3D="50 0 0" Tail3D="0 0 0" ArrowheadHead="Full" AngularSize="120"/></page></CDXML>`;
  const { notes } = fromCdxml(xml);
  const all = notes.join(' ');
  if (!/Generic labels \(R\)/.test(all) || !/hydrogen bonds/.test(all) || !/curved arrows/.test(all)) fail('notes', all);
  else console.log('ok  ', 'notes'.padEnd(20), all);
  for (const broken of ['<CDXML><page>', 'not xml at all', '<?xml version="1.0"?><svg/>']) {
    let threw = false;
    try { fromCdxml(broken); } catch { threw = true; }
    if (!threw) fail('malformed', `accepted ${JSON.stringify(broken)}`);
  }
}

console.log(`${Object.keys(CASES).length + 5} checks, ${bad} failures`);
process.exit(bad ? 1 : 0);
