// File access that works both in the browser and inside the Tauri desktop app.

interface TauriInternals {
  invoke: (cmd: string, args?: Record<string, unknown>) => Promise<unknown>;
}

const tauri = (): TauriInternals | null =>
  (window as unknown as { __TAURI_INTERNALS__?: TauriInternals }).__TAURI_INTERNALS__ ?? null;

export const isDesktop = () => tauri() !== null;
export const isMac = /Mac|iPhone|iPad/.test(navigator.platform) || navigator.userAgent.includes('Mac OS');

export interface FileFilter {
  name: string;
  extensions: string[];
}

/** Save text or bytes. Desktop: native Save panel. Browser: download. Returns false if cancelled. */
export async function saveFile(defaultName: string, data: string | Uint8Array, filter: FileFilter): Promise<boolean> {
  const t = tauri();
  if (t) {
    const path = (await t.invoke('plugin:dialog|save', {
      options: { defaultPath: defaultName, filters: [filter] },
    })) as string | null;
    if (!path) return false;
    const bytes = typeof data === 'string' ? new TextEncoder().encode(data) : data;
    await t.invoke('write_file', { path, contents: Array.from(bytes) });
    return true;
  }
  const blob = new Blob([data as BlobPart], { type: mimeFor(filter.extensions[0]) });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = defaultName;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
  return true;
}

/** Open a web link in the default browser (the desktop app's window never navigates away). */
export function openExternal(url: string) {
  const t = tauri();
  if (t) t.invoke('open_url', { url }).catch(() => {});
  else window.open(url, '_blank', 'noopener,noreferrer');
}

/** Ask the user for a text file. Works with the standard file input in browsers and in the macOS webview. */
export function openTextFile(accept: string): Promise<{ name: string; text: string } | null> {
  return new Promise((resolve) => {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = accept;
    input.onchange = async () => {
      const file = input.files?.[0];
      if (!file) return resolve(null);
      resolve({ name: file.name, text: await file.text() });
    };
    input.click();
  });
}

function mimeFor(ext: string): string {
  switch (ext) {
    case 'png': return 'image/png';
    case 'svg': return 'image/svg+xml';
    case 'cdxml': return 'chemical/x-cdxml';
    case 'json':
    case 'mcv': return 'application/json';
    default: return 'text/plain';
  }
}
