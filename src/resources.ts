// OpenChemLib's force-field and torsion tables (~1.3 MB) are only needed for 3D export, so load them on first use.
import { Resources } from 'openchemlib';
import resourcesUrl from '../node_modules/openchemlib/dist/resources.json?url';

let loading: Promise<void> | null = null;
export function loadChemResources(): Promise<void> {
  loading ??= Resources.registerFromUrl(resourcesUrl).catch((e) => { loading = null; throw e; });
  return loading;
}
