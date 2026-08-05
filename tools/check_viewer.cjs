const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');

const baseUrl = process.argv[2] || 'http://127.0.0.1:8090/viewer/';
const screenshots = path.join(__dirname, '..', 'viewer', 'screenshots');
fs.mkdirSync(screenshots, { recursive: true });
let activeBrowser = null;

async function waitUntilLoaded(page, model) {
  await page.waitForFunction(expected => {
    const loading = document.querySelector('#loading');
    const diagnostics = window.__viewerDiagnostics?.();
    return loading?.classList.contains('hidden') && diagnostics?.model === expected && diagnostics.drawCalls > 0;
  }, model, { timeout: 60000 });
  await page.waitForTimeout(250);
}

(async () => {
  const browser = await chromium.launch({ channel: 'msedge', headless: true });
  activeBrowser = browser;
  const errors = [];
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 }, deviceScaleFactor: 1 });
  page.on('console', message => { if (message.type() === 'error') errors.push(message.text()); });
  page.on('pageerror', error => errors.push(error.message));

  await page.goto(`${baseUrl}?model=1`, { waitUntil: 'domcontentloaded' });
  await waitUntilLoaded(page, 1);
  const split = await page.evaluate(() => window.__viewerDiagnostics());
  const splitPixels = await page.evaluate(() => window.__viewerPixelSample());
  if (split.mode !== 'split' || split.renderTriangles < 140 || splitPixels.uniqueColors < 3)
    throw new Error(`desktop split rendering failed: ${JSON.stringify({ split, splitPixels })}`);

  await page.click('[data-mode="overlay"]');
  await page.waitForFunction(() => window.__viewerDiagnostics().mode === 'overlay');
  const overlay = await page.evaluate(() => window.__viewerDiagnostics());
  const overlayPixels = await page.evaluate(() => window.__viewerPixelSample());
  if (!overlay.originalVisible || !overlay.proxyVisible || overlayPixels.uniqueColors < 3)
    throw new Error(`overlay rendering failed: ${JSON.stringify(overlay)}`);

  await page.selectOption('#modelSelect', '16');
  await waitUntilLoaded(page, 16);
  await page.click('[data-mode="split"]');
  const large = await page.evaluate(() => window.__viewerDiagnostics());
  const largePixels = await page.evaluate(() => window.__viewerPixelSample());
  if (large.renderTriangles < 164597 || large.cameraDistance <= 0 || largePixels.uniqueColors < 3)
    throw new Error(`large model rendering failed: ${JSON.stringify(large)}`);
  await page.screenshot({ path: path.join(screenshots, 'desktop-model-16.png'), fullPage: true });

  const layout = await page.evaluate(() => ({
    overflow: document.documentElement.scrollWidth > innerWidth,
    panel: document.querySelector('#panel').getBoundingClientRect().toJSON(),
    labels: document.querySelector('#labels').getBoundingClientRect().toJSON()
  }));
  if (layout.overflow || layout.panel.right > 1440 || layout.panel.bottom > 900)
    throw new Error(`desktop layout failed: ${JSON.stringify(layout)}`);

  const mobile = await browser.newPage({ viewport: { width: 390, height: 844 }, deviceScaleFactor: 1 });
  mobile.on('console', message => { if (message.type() === 'error') errors.push(message.text()); });
  mobile.on('pageerror', error => errors.push(error.message));
  await mobile.goto(`${baseUrl}?model=6`, { waitUntil: 'domcontentloaded' });
  await waitUntilLoaded(mobile, 6);
  const mobileState = await mobile.evaluate(() => ({
    diagnostics: window.__viewerDiagnostics(),
    pixels: window.__viewerPixelSample(),
    overflow: document.documentElement.scrollWidth > innerWidth,
    panel: document.querySelector('#panel').getBoundingClientRect().toJSON()
  }));
  if (mobileState.overflow || mobileState.panel.left < 0 || mobileState.panel.right > 390 ||
      mobileState.panel.bottom > 844 || mobileState.pixels.uniqueColors < 3)
    throw new Error(`mobile rendering/layout failed: ${JSON.stringify(mobileState)}`);
  await mobile.screenshot({ path: path.join(screenshots, 'mobile-model-6.png'), fullPage: true });

  if (errors.length) throw new Error(`browser errors: ${errors.join(' | ')}`);
  console.log(JSON.stringify({ ok: true, split, overlay, large, layout, mobileState }));
  await browser.close();
  activeBrowser = null;
})().catch(error => {
  console.error(error);
  if (activeBrowser) activeBrowser.close().catch(() => {});
  process.exitCode = 1;
});
