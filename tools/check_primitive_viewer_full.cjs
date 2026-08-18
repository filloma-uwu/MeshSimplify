const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');

const baseUrl = process.argv[2] || 'http://127.0.0.1:8137/viewer/primitive_analysis.html?manifest=/outputs/direct_coverage_uniform_closed_responsibility_full_v13/viewer_manifest.json&model=2';
const screenshots = path.join(__dirname, '..', 'viewer', 'screenshots');
fs.mkdirSync(screenshots, { recursive: true });
let activeBrowser = null;

async function waitLoaded(page, model) {
  await page.waitForFunction(({ expected, strength }) => {
    const state = window.__primitiveViewerDiagnostics?.();
    return state && !state.loading && state.model === expected && state.drawCalls > 0 &&
      (!state.strengthVariantCount || Math.abs(state.analysisStrength - strength) < 1e-9);
  }, { expected: model, strength: 0.41 }, { timeout: 180000 });
  await page.waitForTimeout(250);
}

(async () => {
  const browser = await chromium.launch({ channel: 'msedge', headless: true });
  activeBrowser = browser;
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const errors = [];
  page.on('console', message => { if (message.type() === 'error') errors.push(message.text()); });
  page.on('pageerror', error => errors.push(error.message));
  await page.goto(baseUrl, { waitUntil: 'domcontentloaded' });
  await waitLoaded(page, 2);
  const modelIds = await page.locator('#modelSelect option').evaluateAll(options => options.map(option => Number(option.value)));
  const expected = [2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20];
  if (JSON.stringify(modelIds) !== JSON.stringify(expected)) throw new Error(`manifest models mismatch: ${modelIds}`);

  const results = [];
  for (const modelId of modelIds) {
    await page.selectOption('#modelSelect', String(modelId));
    await waitLoaded(page, modelId);
    const check = await page.evaluate(() => ({
      state: window.__primitiveViewerDiagnostics(),
      pixels: window.__primitiveViewerPixelSample(),
    }));
    if (!check.state.manifestComplete || check.state.manifestModelCount !== 13 ||
        check.state.primitiveMeshes < 1 ||
        check.state.visiblePrimitives !== check.state.primitiveCount ||
        check.state.sourceMeshes < 1 || check.state.regionMeshes < 1 ||
        check.state.regionMeshes > check.state.primitiveCount ||
        check.state.cameraDistance <= 0 || check.pixels.uniqueColors < 3 || check.pixels.nonBackground < 20) {
      throw new Error(`model ${modelId} failed: ${JSON.stringify(check)}`);
    }
    results.push({
      model: modelId,
      primitives: check.state.primitiveCount,
      triangles: check.state.triangles,
      colors: check.pixels.uniqueColors,
    });
  }

  await page.selectOption('#modelSelect', '16');
  await waitLoaded(page, 16);
  await page.click('[data-mode="phase3"]');
  await page.waitForFunction(() => window.__primitiveViewerDiagnostics().visiblePrimitives > 0);
  await page.click('[data-mode="split"]');
  await page.waitForFunction(() => window.__primitiveViewerDiagnostics().mode === 'split');
  await page.screenshot({ path: path.join(screenshots, 'primitive-analysis-full-model-16.png'), fullPage: true });

  await page.selectOption('#modelSelect', '18');
  await waitLoaded(page, 18);
  await page.click('[data-mode="error"]');
  await page.waitForFunction(() => window.__primitiveViewerDiagnostics().mode === 'error');
  const errorOverlay = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  if (errorOverlay.maximumErrorDistance > 0 &&
      (errorOverlay.maximumErrorSegmentCount !== 1 ||
       errorOverlay.maximumErrorEndpointCount !== 2 ||
       errorOverlay.maximumErrorLinePixels < 6 ||
       !errorOverlay.maximumErrorOverlayDepthTest ||
       errorOverlay.maximumErrorOverlayOrder < 100)) {
    throw new Error(`maximum-error overlay failed: ${JSON.stringify(errorOverlay)}`);
  }
  await page.screenshot({ path: path.join(screenshots, 'primitive-analysis-error-model-18.png'), fullPage: true });

  const mobile = await browser.newPage({ viewport: { width: 390, height: 844 } });
  mobile.on('console', message => { if (message.type() === 'error') errors.push(message.text()); });
  mobile.on('pageerror', error => errors.push(error.message));
  const mobileUrl = new URL(baseUrl);
  mobileUrl.searchParams.set('model', '5');
  await mobile.goto(mobileUrl.toString(), { waitUntil: 'domcontentloaded' });
  await waitLoaded(mobile, 5);
  const mobileState = await mobile.evaluate(() => ({
    diagnostics: window.__primitiveViewerDiagnostics(),
    pixels: window.__primitiveViewerPixelSample(),
    overflow: document.documentElement.scrollWidth > innerWidth,
    panel: document.querySelector('#panel').getBoundingClientRect().toJSON(),
  }));
  if (mobileState.overflow || mobileState.panel.left < 0 || mobileState.panel.right > 390 ||
      mobileState.panel.bottom > 844 || mobileState.pixels.nonBackground < 10) {
    throw new Error(`mobile failed: ${JSON.stringify(mobileState)}`);
  }
  await mobile.screenshot({ path: path.join(screenshots, 'primitive-analysis-full-mobile.png'), fullPage: true });
  if (errors.length) throw new Error(`browser errors: ${errors.join(' | ')}`);
  console.log(JSON.stringify({ ok: true, results, mobile: mobileState.diagnostics }));
  await browser.close();
  activeBrowser = null;
})().catch(error => {
  console.error(error);
  if (activeBrowser) activeBrowser.close().catch(() => {});
  process.exitCode = 1;
});
