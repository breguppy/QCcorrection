const { app, BrowserWindow } = require('electron');
const path = require('path');
const { spawn } = require('child_process');
const http = require('http');

const PORT = 3170; // fixed localhost port
let rProc, win;

function rscriptPath() {
  if (process.platform === 'win32') {
    return path.join(process.resourcesPath, 'r-env', 'win', 'R', 'bin', 'Rscript.exe');
  } else {
    return path.join(process.resourcesPath, 'r-env', 'mac', 'R.framework', 'Resources', 'bin', 'Rscript');
  }
}

function envForR() {
  const res = process.resourcesPath;
  const isWin = process.platform === 'win32';
  const libDir = path.join(res, 'r-env', isWin ? 'win' : 'mac', 'library');
  const pandoc = path.join(res, 'r-env', isWin ? 'win' : 'mac', 'pandoc');

  const env = Object.assign({}, process.env, {
    SHINY_PORT: String(PORT),
    R_PACK_LIB: libDir,
    PANDOC_PATH: pandoc
  });
  if (isWin) {
    env.PATH = path.join(res, 'r-env', 'win', 'R', 'bin', 'x64') + ';' + env.PATH;
  }
  return env;
}

function waitForServer(url, timeoutMs = 20000, intervalMs = 200) {
  const start = Date.now();
  return new Promise((resolve, reject) => {
    (function tick() {
      http.get(url, res => {
        res.resume();
        if (res.statusCode >= 200 && res.statusCode < 500) return resolve();
        setTimeout(tick, intervalMs);
      }).on('error', () => setTimeout(tick, intervalMs));
      if (Date.now() - start > timeoutMs) reject(new Error('Shiny did not start in time'));
    })();
  });
}

async function createWindow() {
  const startR = path.join(process.resourcesPath, 'scripts', 'start.R');

  rProc = spawn(rscriptPath(), [startR], { env: envForR(), windowsHide: true });
  rProc.on('exit', () => { if (win && !win.isDestroyed()) win.close(); });

  await waitForServer(`http://127.0.0.1:${PORT}/`);

  win = new BrowserWindow({ width: 1280, height: 800, webPreferences: { contextIsolation: true, nodeIntegration: false } });
  win.removeMenu();
  await win.loadURL(`http://127.0.0.1:${PORT}/`);
}

app.whenReady().then(createWindow);
app.on('window-all-closed', () => { if (process.platform !== 'darwin') app.quit(); });
app.on('before-quit', () => { try { rProc?.kill(); } catch {} });
