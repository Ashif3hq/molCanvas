import { defineConfig } from 'vite';

// base './' keeps asset paths relative so the same build works on GitHub Pages and inside Tauri.
export default defineConfig({
  base: './',
  clearScreen: false,
  server: { port: 1420, strictPort: true },
  build: { target: 'safari15', chunkSizeWarningLimit: 4000 },
});
