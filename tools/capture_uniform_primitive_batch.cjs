const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');

const baseUrl = process.argv[2] ||
  'http://127.0.0.1:8093/viewer/primitive_analysis.html' +
  '?manifest=/outputs/surface_primitives_v190_uniform_final/viewer_manifest.json';
const outputDir = path.resolve(process.argv[3] ||
  path.join(__dirname, '..', 'viewer', 'screenshots', 'v190'));
const requestedMode = process.argv[4] || 'split';
fs.mkdirSync(outputDir, { recursive: true });

async function waitLoaded(page, model) {
  await page.waitForFunction(expected => {
    const state = window.__primitiveViewerDiagnostics?.();
    return state?.model === expected && state?.manifestComplete &&
      !state?.loading;
  }, model, { timeout: 120000 });
  await page.waitForTimeout(250);
}

let browser;
(async () => {
  browser = await chromium.launch({ headless: true });
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const errors = [];
  page.on('console', message => {
    if (message.type() === 'error') errors.push(message.text());
  });
  page.on('pageerror', error => errors.push(error.message));

  await page.goto(baseUrl, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => {
    const state = window.__primitiveViewerDiagnostics?.();
    return state?.manifestComplete && state.model !== null;
  });
  const first = await page.evaluate(() => window.__primitiveViewerDiagnostics().model);
  await waitLoaded(page, first);
  const modelIds = await page.locator('#modelSelect option').evaluateAll(
    options => options.map(option => Number(option.value)));

  const records = [];
  for (const modelId of modelIds) {
    console.log(`checking_model=${modelId}`);
    if (modelId !== first) await page.selectOption('#modelSelect', String(modelId));
    await waitLoaded(page, modelId);
    await page.click(`[data-mode="${requestedMode}"]`);
    await page.click('#fitButton');

    let state = await page.evaluate(() => window.__primitiveViewerDiagnostics());
    if (state.primitiveOptionCount !== state.primitiveCount ||
        state.selectedPrimitive !== 'all' || !state.allBaseMaterialsOpaque) {
      throw new Error(`initial state failed for model ${modelId}: ${JSON.stringify(state)}`);
    }
    const pixels = await page.evaluate(() => window.__primitiveViewerPixelSample());
    if (pixels.uniqueColors < 3 || pixels.nonBackground < 20) {
      throw new Error(`blank canvas for model ${modelId}: ${JSON.stringify(pixels)}`);
    }

    const last = String(state.primitiveCount - 1);
    await page.selectOption('#primitiveSelect', last);
    await page.waitForFunction(expected =>
      window.__primitiveViewerDiagnostics?.().selectedPrimitive === expected, last);
    state = await page.evaluate(() => window.__primitiveViewerDiagnostics());
    if (state.selectionOverlayCount === 0 && state.selectedBaseMeshCount === 0) {
      throw new Error(`primitive ${last} is not selectable for model ${modelId}`);
    }
    await page.selectOption('#primitiveSelect', 'all');
    await page.waitForFunction(() => {
      const current = window.__primitiveViewerDiagnostics?.();
      return current?.selectedPrimitive === 'all' && current.allBaseMaterialsOpaque;
    });
    await page.screenshot({
      path: path.join(outputDir, `model-${modelId}.png`),
      fullPage: true,
    });
    records.push({ modelId, primitiveCount: state.primitiveCount, pixels });
  }

  if (errors.length) throw new Error(`browser errors: ${errors.join(' | ')}`);
  console.log(JSON.stringify({ ok: true, models: records }, null, 2));
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
}).finally(async () => {
  if (browser) await browser.close();
});
