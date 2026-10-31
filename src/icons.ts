// Inline SVG icons (24×24, stroke = currentColor).
const svg = (body: string, extra = '') =>
  `<svg viewBox="0 0 24 24" width="22" height="22" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"${extra}>${body}</svg>`;

const ring = (n: number, aromatic = false) => {
  const r = 8, cx = 12, cy = 12.5;
  const pts: [number, number][] = [];
  for (let i = 0; i < n; i++) {
    const a = -Math.PI / 2 + (2 * Math.PI * i) / n;
    pts.push([cx + r * Math.cos(a), cy + r * Math.sin(a)]);
  }
  let body = `<polygon points="${pts.map((p) => p.map((v) => v.toFixed(2)).join(',')).join(' ')}"/>`;
  if (aromatic) {
    for (let i = 0; i < n; i += 2) {
      const [x1, y1] = pts[i], [x2, y2] = pts[(i + 1) % n];
      const k = 0.72;
      body += `<line x1="${(cx + (x1 - cx) * k).toFixed(2)}" y1="${(cy + (y1 - cy) * k).toFixed(2)}" x2="${(cx + (x2 - cx) * k).toFixed(2)}" y2="${(cy + (y2 - cy) * k).toFixed(2)}"/>`;
    }
  }
  return svg(body);
};

export const ICONS: Record<string, string> = {
  select: svg('<path d="M5 3l13 8-6 1.5L9 19z" fill="currentColor" fill-opacity="0.15"/>'),
  erase: svg('<path d="M15.5 4.5l4 4L10 18H6l-2-2z"/><path d="M11 9l4 4"/><path d="M10 20h10"/>'),
  single: svg('<line x1="5" y1="17" x2="19" y2="7"/>'),
  double: svg('<line x1="4.5" y1="15" x2="17.5" y2="5.5"/><line x1="6.5" y1="18.5" x2="19.5" y2="9"/>'),
  triple: svg('<line x1="4" y1="14" x2="16.5" y2="5"/><line x1="5.5" y1="17" x2="18" y2="8"/><line x1="7" y1="20" x2="19.5" y2="11"/>'),
  wedge: svg('<polygon points="5,17 18,5.5 20,8.5" fill="currentColor"/>'),
  hash: svg('<g stroke-width="1.5"><line x1="5.5" y1="16.2" x2="5.9" y2="16.8"/><line x1="8.3" y1="13.5" x2="9.3" y2="15"/><line x1="11" y1="10.9" x2="12.6" y2="13.2"/><line x1="13.8" y1="8.3" x2="15.9" y2="11.4"/><line x1="16.6" y1="5.7" x2="19.2" y2="9.6"/></g>'),
  wavy: svg('<path d="M4 17c1.5-3 3.5 0 5-3s3.5 0 5-3 3.5 0 5-3"/>'),
  chain: svg('<polyline points="3,15 7.5,9 12,15 16.5,9 21,15"/>'),
  ring3: ring(3), ring4: ring(4), ring5: ring(5), ring6: ring(6), ring7: ring(7), ring8: ring(8),
  benzene: ring(6, true),
  plus: svg('<circle cx="12" cy="12" r="8"/><path d="M12 8v8M8 12h8"/>'),
  minus: svg('<circle cx="12" cy="12" r="8"/><path d="M8 12h8"/>'),
  radical: svg('<line x1="4.5" y1="19.5" x2="13" y2="11"/><circle cx="17.5" cy="6.5" r="2.4" fill="currentColor" stroke="none"/>'),
  lonepairs: svg('<circle cx="12" cy="13.5" r="4.2"/><g fill="currentColor" stroke="none"><circle cx="10" cy="5.6" r="1.3"/><circle cx="14" cy="5.6" r="1.3"/><circle cx="4.4" cy="11.5" r="1.3"/><circle cx="4.4" cy="15.5" r="1.3"/><circle cx="19.6" cy="11.5" r="1.3"/><circle cx="19.6" cy="15.5" r="1.3"/></g>'),
  arrow: svg('<line x1="3" y1="12" x2="19" y2="12"/><polygon points="21,12 16,9 16,15" fill="currentColor"/>'),
  equilibrium: svg('<line x1="3" y1="9.5" x2="20" y2="9.5"/><path d="M20 9.5l-4-3"/><line x1="4" y1="14.5" x2="21" y2="14.5"/><path d="M4 14.5l4 3"/>'),
  retro: svg('<line x1="3" y1="9.5" x2="17" y2="9.5"/><line x1="3" y1="14.5" x2="17" y2="14.5"/><path d="M14 5l7 7-7 7"/>'),
  resonance: svg('<line x1="6" y1="12" x2="18" y2="12"/><polygon points="21,12 16,9 16,15" fill="currentColor"/><polygon points="3,12 8,9 8,15" fill="currentColor"/>'),
  text: svg('<path d="M5 6V4h14v2M12 4v16M9 20h6"/>'),
  undo: svg('<path d="M9 14L4 9l5-5"/><path d="M4 9h10a6 6 0 010 12h-3"/>'),
  redo: svg('<path d="M15 14l5-5-5-5"/><path d="M20 9H10a6 6 0 000 12h3"/>'),
  new: svg('<path d="M14 3H6a1 1 0 00-1 1v16a1 1 0 001 1h12a1 1 0 001-1V8z"/><path d="M14 3v5h5M12 11v6M9 14h6"/>'),
  open: svg('<path d="M3 7a2 2 0 012-2h4l2 2h8a2 2 0 012 2v8a2 2 0 01-2 2H5a2 2 0 01-2-2z"/>'),
  save: svg('<path d="M5 3h11l5 5v11a2 2 0 01-2 2H5a2 2 0 01-2-2V5a2 2 0 012-2z"/><path d="M7 3v5h8V3M7 21v-7h10v7"/>'),
  export: svg('<path d="M12 15V3M7 8l5-5 5 5"/><path d="M5 13v6a2 2 0 002 2h10a2 2 0 002-2v-6"/>'),
  copy: svg('<rect x="8" y="8" width="12" height="12" rx="2"/><path d="M16 8V6a2 2 0 00-2-2H6a2 2 0 00-2 2v8a2 2 0 002 2h2"/>'),
  clean: svg('<path d="M12 3l1.8 4.7L18.5 9l-4.7 1.8L12 15.5l-1.8-4.7L5.5 9l4.7-1.3z"/><path d="M18 15l.8 2 2 .8-2 .8-.8 2-.8-2-2-.8 2-.8z"/>'),
  rotate: svg('<path d="M20 12a8 8 0 11-2.3-5.6"/><path d="M20 4v5h-5"/>'),
  flipH: svg('<path d="M12 3v18" stroke-dasharray="2 2"/><path d="M9 7L4 17h5zM15 7l5 10h-5z"/>'),
  flipV: svg('<path d="M3 12h18" stroke-dasharray="2 2"/><path d="M7 9L17 4v5zM7 15l10 5v-5z"/>'),
  zoomIn: svg('<circle cx="11" cy="11" r="7"/><path d="M21 21l-5-5M11 8v6M8 11h6"/>'),
  zoomOut: svg('<circle cx="11" cy="11" r="7"/><path d="M21 21l-5-5M8 11h6"/>'),
  fit: svg('<path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5"/>'),
  sun: svg('<circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/>'),
  moon: svg('<path d="M20 14.5A8 8 0 019.5 4a8 8 0 1010.5 10.5z"/>'),
  palette: svg('<path d="M12 3a9 9 0 100 18c1.2 0 1.5-.8 1.5-1.5 0-.9-.7-1.2-.7-2 0-.8.6-1.5 1.5-1.5H17a4 4 0 004-4c0-5-4-9-9-9z"/><circle cx="7.5" cy="11" r="1.2" fill="currentColor"/><circle cx="10.5" cy="7" r="1.2" fill="currentColor"/><circle cx="15" cy="7.5" r="1.2" fill="currentColor"/>'),
  search: svg('<circle cx="11" cy="11" r="7"/><path d="M21 21l-5-5"/>'),
  library: svg('<rect x="4" y="4" width="6" height="6" rx="1"/><rect x="14" y="4" width="6" height="6" rx="1"/><rect x="4" y="14" width="6" height="6" rx="1"/><rect x="14" y="14" width="6" height="6" rx="1"/>'),
  table: svg('<rect x="3" y="5" width="18" height="14" rx="2"/><path d="M3 10h18M3 14.5h18M8 5v14M13 5v14M18 5v14"/>'),
  bracket: svg('<path d="M7 4H4v16h3M13 4h3v16h-3"/><text x="17.5" y="22" font-size="6.5" fill="currentColor" stroke="none" font-style="italic">m</text>'),
  polymer: svg('<path d="M6 5c-2 3-2 11 0 14M18 5c2 3 2 11 0 14"/><polyline points="8,14 10.5,10 13.5,14 16,10"/><text x="19" y="22" font-size="6" fill="currentColor" stroke="none" font-style="italic">n</text>'),
};
