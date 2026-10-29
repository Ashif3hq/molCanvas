MolCanvas for Mac: how to install
=================================

MolCanvas is a free chemical structure editor for organic and polymer chemists, by Ashif.
This build runs natively on Apple Silicon Macs (M1, M2, M3, M4, M5) with macOS 11 or newer.

1. Download
   Download the file ending in "_aarch64.dmg" from the Assets list on this release page.

2. Install
   Double-click the .dmg file. In the window that opens, drag MolCanvas onto the
   Applications folder. Then eject the MolCanvas disk (the eject arrow next to it in Finder).

3. Open it the first time
   MolCanvas is not notarised by Apple yet, so macOS asks you to confirm the first launch:
   - Open the Applications folder, right-click (or Control-click) MolCanvas and choose Open,
     then click Open in the dialog.
   - On macOS 15 (Sequoia) or newer, if you only see "Done" or "Move to Bin": click Done,
     open System Settings > Privacy & Security, scroll down and click "Open Anyway" next to
     the MolCanvas message, then confirm with your password or Touch ID.
   - If macOS says the app "is damaged and can't be opened", open Terminal and run:
         xattr -dr com.apple.quarantine /Applications/MolCanvas.app
     then open MolCanvas again. Nothing is wrong with the app: this removes the
     "downloaded from the internet" flag that macOS puts on apps that are not notarised.
   After the first launch, MolCanvas opens normally like any other app.

4. Getting started
   - Pick a bond, ring or atom tool on the left and click on the canvas to draw.
   - Hover an atom and press a letter to change the element (C, N, O, S, F, L for Cl,
     B for Br), + or - for a charge, or . for a radical.
   - The library on the left has 350+ ready-made molecules, monomers and polymers.
   - The bottom info bar explains every atom and bond, and why a change is refused
     when it breaks valence rules.
   - Open and save .mcv, ChemDraw .cdxml, .mol, .sdf and SMILES. Export PNG, SVG, 3D SDF,
     MOL2, PDB and XYZ.

Updating: download the newer .dmg and drag MolCanvas to Applications again, choosing Replace.
Uninstalling: drag MolCanvas from Applications to the Bin.

Try it without installing: https://ashif3hq.github.io/molCanvas/
Source code and issues: https://github.com/Ashif3hq/molCanvas
