const { app, BrowserWindow } = require('electron');
const path = require('path');
const fs = require('fs');
const { spawn } = require('child_process');
const http = require('http');

const PORT = 3170;
let rProc, win;
const isDev = !app.isPackaged;
const BASE = isDev ? path.join(__dirname, '..') : process.resourcesPath;

function findRscript() {
  const cands = [
    path.join(BASE, 'r-env', 'win', 'R', 'bin', 'x64', 'Rscript.exe'),
    path.join(BASE, 'r-env', 'win', 'R', 'R-4.4.1', 'bin', 'x64', 'Rscript.exe'),
    path.join(BASE, 'r-env', 'win', 'R', 'R-4.4.2', 'bin', 'x64', 'Rscript.exe')
  ];
  for (const p of cands) if (fs.existsSync(p)) return p;
  throw new Error('Rscript.exe not found. Checked:\n' + cands.join('\n'));
}

function envForR() {
  const base = BASE;
  const libDir = path.join(base, 'r-env', 'win', 'library');

  // find pandoc.exe
  const p1 = path.join(base, 'r-env', 'win', 'pandoc');                  // unversioned
  const p2 = path.join(base, 'r-env', 'win', 'pandoc', 'pandoc-3.8.2');   // versioned
  const pandocDir = fs.existsSync(path.join(p1, 'pandoc.exe')) ? p1 :
                    fs.existsSync(path.join(p2, 'pandoc.exe')) ? p2 : '';

  const rscript = findRscript();
  const rBin = path.dirname(rscript);

  return {
    ...process.env,
    SHINY_PORT: String(PORT),
    R_PACK_LIB: libDir,
    RSTUDIO_PANDOC: pandocDir,
    ELECTRON: "1",
    PATH: `${rBin};${process.env.PATH || ''}`
  };
}


function waitForServer(url, timeoutMs = 60000, intervalMs = 200) {
  const t0 = Date.now();
  return new Promise((resolve, reject) => {
    (function poll() {
      http.get(url, res => {
        res.resume();
        if (res.statusCode >= 200 && res.statusCode < 500) return resolve();
        setTimeout(poll, intervalMs);
      }).on('error', () => setTimeout(poll, intervalMs));
      if (Date.now() - t0 > timeoutMs) reject(new Error('Shiny did not start in time'));
    })();
  });
}

async function createWindow() {
  const rscript = findRscript();
  const env = envForR(rscript);
  const startR = path.join(BASE, 'scripts', 'start.R');

  rProc = spawn(rscript, [startR], { env, windowsHide: true });
  rProc.stdout?.on('data', d => console.log('[R]', d.toString()));
  rProc.stderr?.on('data', d => console.error('[R]', d.toString()));
  rProc.on('exit', c => console.log('R exited with', c));

  await waitForServer(`http://127.0.0.1:${PORT}/`);
  win = new BrowserWindow({ width: 1280, height: 800, webPreferences: { contextIsolation: true, nodeIntegration: false } });
  win.removeMenu();
  await win.loadURL(`http://127.0.0.1:${PORT}/`);
}

app.whenReady().then(createWindow);
app.on('window-all-closed', () => { if (process.platform !== 'darwin') app.quit(); });
app.on('before-quit', () => { try { rProc?.kill(); } catch {} });
