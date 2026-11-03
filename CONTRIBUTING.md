# Contributing to MolCanvas

Thanks for your interest in MolCanvas, a free chemical structure editor for organic and polymer
chemists, created and maintained by **Ashif**.

Please read this page before opening an issue or a pull request. Everyone taking part is expected
to follow the [Code of Conduct](CODE_OF_CONDUCT.md).

## Before you start: the licence

MolCanvas is **free to use but all rights reserved** (see [LICENSE](LICENSE)). It is not an open
source project in the usual sense, which affects contributions:

- Bug reports, feature ideas, preset requests and feedback in issues are always welcome.
- Code, artwork and documentation changes are accepted only through pull requests to this
  repository, and only ones Ashif has agreed to. **Open an issue first** and wait for a go-ahead
  before writing anything larger than a small fix, so your time isn't wasted.
- Forking is allowed only as far as GitHub needs to make a pull request. A fork may not be
  published, redistributed or hosted elsewhere (see LICENSE).

**By opening a pull request you confirm that the work is your own (or you have the right to
submit it), and you grant Ashif a perpetual, worldwide, royalty-free, irrevocable licence to use,
modify and distribute it as part of MolCanvas, under MolCanvas's own licence terms.** This grant
cannot be used to relicense your contribution, or MolCanvas, under any other licence. You will be
credited in the commit history and the release notes.

## Reporting bugs

Search the [existing issues](https://github.com/Ashif3hq/molCanvas/issues) first. A good report
includes:

- the MolCanvas version, and whether you used the Mac app or the
  [web version](https://ashif3hq.github.io/molCanvas/) (plus your macOS or browser version);
- what you did, what you expected, and what happened instead;
- the structure involved, as SMILES or an attached `.mol`, `.sdf` or `.cdxml` file, and a
  screenshot if the problem is visual.

**Chemistry mistakes count as bugs.** If MolCanvas allows an impossible structure, gives a wrong
formula, mass, name or solubility value, or writes a file another program reads differently,
please say what the correct result is and where it comes from (a textbook, paper, database entry
or another program).

**Security problems must not go in public issues.** Follow [SECURITY.md](SECURITY.md) instead.

## Suggesting features and presets

Open an issue describing the problem it solves and who it helps (for example, polymer chemists
drawing block copolymers). For a new preset molecule, monomer or repeat unit, give its name and a
SMILES string, ideally with a PubChem or literature reference. Check the
[roadmap](README.md#roadmap) first; it may already be planned.

## Working on the code

Requirements: Node.js 20+, and for the desktop app [Rust](https://rustup.rs) plus the Xcode
command line tools on a Mac with Apple Silicon.

```sh
npm install
npm run dev        # web version at http://localhost:1420
npm run app:dev    # desktop app with hot reload (macOS)
```

The [README](README.md#how-it-works) explains what each file in `src/` does.

### Ground rules

- **Chemistry must be correct.** MolCanvas follows chemical theory: valence rules, charges and
  radicals are enforced, and anything it can't represent is explained in the info bar instead of
  being changed silently. New features must keep it that way. When a calculation is not
  supported for a structure, refuse and say why rather than return a guess.
- **Keep it small and focused.** One fix or feature per pull request, written in TypeScript in the
  style of the surrounding code. Avoid new dependencies unless agreed in the issue.
- **Explain things to the user.** User-facing text should be plain English that a chemistry
  student understands.
- **Both versions must work.** Changes should work in the web version and in the Mac app.

### Checks to run before opening a pull request

```sh
npm run build              # type-checks and builds the web version
npm run check:presets      # every preset structure is valid
npm run check:solubility   # group-contribution sums against hand-worked values
npm run check:hansen       # Hansen distance, mixtures and sphere fitting
npm run check:export       # MOL and 3D export round trips
npm run check:cdxml        # ChemDraw CDXML round trips
npm run check:radicals     # radicals, carbenes and ions
npm run check:brackets     # polymer brackets
```

If you change a calculation or a file format, add or update a check in `scripts/` that proves the
new behaviour, with reference values and where they come from.

### Pull requests

1. Link the issue the pull request addresses.
2. Describe what changed and how you tested it; add before/after screenshots for visual changes.
3. Make sure the **macOS app** check on the pull request passes.
4. Don't change version numbers or release files; Ashif handles releases.

Thank you for helping make MolCanvas better.
