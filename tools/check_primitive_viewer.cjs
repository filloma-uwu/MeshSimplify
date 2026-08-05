const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');

const baseUrl = process.argv[2] || 'http://127.0.0.1:8091/viewer/primitive_analysis.html?manifest=/outputs/primitive_analysis_viewer_demo_v2/viewer_manifest.json';
const screenshots = path.join(__dirname, '..', 'viewer', 'screenshots');
fs.mkdirSync(screenshots, { recursive: true });
let activeBrowser = null;

async function waitLoaded(page, model) {
  await page.waitForFunction(expected => {
    const state = window.__primitiveViewerDiagnostics?.();
    return state && !state.loading && state.model === expected && state.drawCalls > 0;
  }, model, { timeout: 60000 });
  await page.waitForTimeout(300);
}

(async () => {
  const browser = await chromium.launch({ channel: 'msedge', headless: true });
  activeBrowser = browser;
  const errors = [];
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  page.on('console', message => {
    if (message.type() === 'error') {
      errors.push(message.text());
      console.error(`browser console: ${message.text()}`);
    }
  });
  page.on('pageerror', error => {
    errors.push(error.message);
    console.error(`browser pageerror: ${error.message}`);
  });
  await page.goto(`${baseUrl}&model=1`, { waitUntil: 'domcontentloaded' });
  await waitLoaded(page, 1);

  const split = await page.evaluate(() => ({ state: window.__primitiveViewerDiagnostics(), pixels: window.__primitiveViewerPixelSample() }));
  if (split.state.mode !== 'split' || split.state.primitiveMeshes !== 14 || !split.state.sourceVisible ||
      !split.state.primitivesVisible || split.pixels.uniqueColors < 3 || split.pixels.nonBackground < 20) {
    throw new Error(`split rendering failed: ${JSON.stringify(split)}`);
  }

  await page.click('[data-mode="regions"]');
  await page.waitForFunction(() => window.__primitiveViewerDiagnostics().mode === 'regions');
  const regions = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  if (regions.visibleRegions !== 14 || regions.primitivesVisible) throw new Error(`regions mode failed: ${JSON.stringify(regions)}`);

  await page.click('[data-mode="primitives"]');
  await page.locator('[data-type="rss"]').uncheck();
  const filtered = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  if (filtered.visiblePrimitives !== 0) throw new Error(`type filter failed: ${JSON.stringify(filtered)}`);
  await page.locator('[data-type="rss"]').check();
  await page.selectOption('#primitiveSelect', '0');
  const selected = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  if (selected.visiblePrimitives !== 1 || selected.selectedPrimitive !== '0') throw new Error(`primitive selection failed: ${JSON.stringify(selected)}`);

  await page.selectOption('#modelSelect', '6');
  await waitLoaded(page, 6);
  await page.locator('[data-type="rss"]').uncheck();
  const capsules = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  if (capsules.primitiveCount !== 19 || capsules.visiblePrimitives !== 2) {
    throw new Error(`real-model capsule filter failed: ${JSON.stringify(capsules)}`);
  }
  await page.locator('[data-type="rss"]').check();
  await page.fill('#opacity', '40');
  await page.check('#wireframe');
  await page.click('#fitButton');

  await page.selectOption('#modelSelect', '21');
  await waitLoaded(page, 21);
  const switched = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  if (switched.primitiveMeshes !== 1 || switched.primitiveCount !== 1 || switched.cameraDistance <= 0 ||
      Math.abs(switched.primitiveOpacity - 0.4) > 1e-6 || !switched.primitiveWireframe) {
    throw new Error(`model switching failed: ${JSON.stringify(switched)}`);
  }
  await page.screenshot({ path: path.join(screenshots, 'primitive-analysis-desktop.png'), fullPage: true });
  const desktopLayout = await page.evaluate(() => ({
    overflow: document.documentElement.scrollWidth > innerWidth,
    panel: document.querySelector('#panel').getBoundingClientRect().toJSON(),
    topbar: document.querySelector('#topbar').getBoundingClientRect().toJSON(),
  }));
  if (desktopLayout.overflow || desktopLayout.panel.right > 1440 || desktopLayout.panel.bottom > 900) {
    throw new Error(`desktop layout failed: ${JSON.stringify(desktopLayout)}`);
  }

  const mobile = await browser.newPage({ viewport: { width: 390, height: 844 } });
  mobile.on('console', message => {
    if (message.type() === 'error') {
      errors.push(message.text());
      console.error(`mobile console: ${message.text()}`);
    }
  });
  mobile.on('pageerror', error => {
    errors.push(error.message);
    console.error(`mobile pageerror: ${error.message}`);
  });
  await mobile.goto(`${baseUrl}&model=1`, { waitUntil: 'domcontentloaded' });
  await waitLoaded(mobile, 1);
  const mobileCheck = await mobile.evaluate(() => ({
    state: window.__primitiveViewerDiagnostics(),
    pixels: window.__primitiveViewerPixelSample(),
    overflow: document.documentElement.scrollWidth > innerWidth,
    panel: document.querySelector('#panel').getBoundingClientRect().toJSON(),
    topbar: document.querySelector('#topbar').getBoundingClientRect().toJSON(),
  }));
  if (mobileCheck.overflow || mobileCheck.panel.left < 0 || mobileCheck.panel.right > 390 ||
      mobileCheck.panel.bottom > 844 || mobileCheck.topbar.right > 390 ||
      mobileCheck.pixels.uniqueColors < 3 || mobileCheck.pixels.nonBackground < 10) {
    throw new Error(`mobile rendering/layout failed: ${JSON.stringify(mobileCheck)}`);
  }
  await mobile.screenshot({ path: path.join(screenshots, 'primitive-analysis-mobile.png'), fullPage: true });

  if (errors.length) throw new Error(`browser errors: ${errors.join(' | ')}`);
  console.log(JSON.stringify({ ok: true, split, regions, filtered, selected, capsules, switched, desktopLayout, mobileCheck }));
  await browser.close();
  activeBrowser = null;
})().catch(error => {
  console.error(error);
  if (activeBrowser) activeBrowser.close().catch(() => {});
  process.exitCode = 1;
});
