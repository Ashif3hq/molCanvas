# Security policy

MolCanvas is made by **Ashif**. Thank you for helping keep it and its users safe.

## Supported versions

Only the latest release gets security fixes. Older versions are not patched; please
update to the newest release instead.

| Version                     | Supported          |
| --------------------------- | ------------------ |
| 0.4.x (latest release)      | :white_check_mark: |
| 0.3.x and earlier           | :x:                |

This covers both ways of using MolCanvas:

- the **Mac app** (Apple Silicon `.dmg`) from this repository's
  [Releases page](https://github.com/Ashif3hq/molCanvas/releases/latest), and
- the **web version** at <https://ashif3hq.github.io/molCanvas/>, which always runs the
  latest code from `main`.

Copies built or hosted by anyone else are not supported.

## Reporting a vulnerability

**Please do not report security problems in public issues, discussions or pull requests.**

Report them privately through GitHub instead:

1. Go to the repository's [Security tab](https://github.com/Ashif3hq/molCanvas/security).
2. Click **Report a vulnerability**
   (or open <https://github.com/Ashif3hq/molCanvas/security/advisories/new> directly).
3. Describe the problem, including:
   - the MolCanvas version, and whether you used the Mac app or the web version
     (plus your macOS or browser version);
   - steps to reproduce it, ideally with a sample file (`.mol`, `.sdf`, `.cdxml`, etc.)
     or the exact input that triggers it;
   - what an attacker could do with it.

Only you and the maintainer can see the report. If you cannot use GitHub's private
reporting, open a public issue that says only that you have a security report and asks for a
private way to send it, without any details of the problem.

## What to expect

- An acknowledgement within **7 days**.
- An assessment and, if the report is confirmed, a plan for a fix, usually within
  **30 days**. MolCanvas is maintained by one person, so complex fixes may take longer; you
  will be kept updated.
- The fix ships in a new release (and on the web version), followed by a GitHub security
  advisory. You will be credited in the advisory unless you would rather not be.

Please give a reasonable chance to release a fix before disclosing the problem publicly.

## Scope

In scope, for example:

- opening or importing a crafted file (MOL, SDF, CDXML, SMILES and other formats) that runs
  code, reads or writes files it should not, or escapes the app's sandbox;
- script injection (XSS) in the web version or the app's window;
- weaknesses in the Mac app's Tauri configuration, such as its permissions or content
  security policy;
- problems in how release builds are produced in this repository's GitHub Actions workflows.

Out of scope:

- bugs that are not security issues, such as wrong chemistry results or crashes with no
  security impact; please open a normal [issue](https://github.com/Ashif3hq/molCanvas/issues)
  for those;
- macOS Gatekeeper warnings about the app not being notarised (this is known and documented
  in the README);
- vulnerabilities in third-party components (OpenChemLib, Tauri, PubChem and their
  dependencies) that do not affect MolCanvas; please report those to their own projects;
- attacks that need physical access to an unlocked Mac, or a device that is already
  compromised.

Good-faith security research that follows this policy is welcome. Testing must not harm
other people, their data, or the services MolCanvas relies on (such as PubChem or GitHub).

## Licence

Reporting a vulnerability does not grant any licence to the MolCanvas source code. See
[LICENSE](LICENSE): MolCanvas is free to use, all rights reserved.
