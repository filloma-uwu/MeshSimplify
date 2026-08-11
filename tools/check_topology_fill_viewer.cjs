const { chromium } = require('playwright');

const url = process.argv[2];
if (!url) throw new Error('usage: node check_topology_fill_viewer.cjs <viewer-url>');

(async () => {
  const browser = await chromium.launch({ headless: true });
  try {
    const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    page.on('console', message => {
      if (message.type() === 'error') errors.push(message.text());
    });
    await page.goto(url, { waitUntil: 'domcontentloaded' });
    await page.waitForFunction(() => {
      const state = window.__primitiveViewerDiagnostics?.();
      return state?.manifestComplete && !state?.loading;
    }, null, { timeout: 120000 });

    const modelIds = await page.locator('#modelSelect option').evaluateAll(
      options => options.map(option => Number(option.value)));
    const results = [];
    for (const modelId of modelIds) {
      await page.selectOption('#modelSelect', String(modelId));
      await page.waitForFunction(expected => {
        const state = window.__primitiveViewerDiagnostics?.();
        return state && !state.loading && state.model === expected && state.drawCalls > 0;
      }, modelId, { timeout: 120000 });

      const stageAudit = await page.evaluate(() => ({
        diagnostics: window.__primitiveViewerDiagnostics(),
        disabled: Object.fromEntries([...document.querySelectorAll('[data-mode]')]
          .map(button => [button.dataset.mode, button.disabled])),
      }));
      const diagnostics = stageAudit.diagnostics;
      if (diagnostics.sourceMeshes < 1 || diagnostics.phase1Meshes < 1 ||
          diagnostics.phase2Meshes !== 0 || diagnostics.primitiveMeshes !== 0 ||
          diagnostics.triangulatedMeshes !== 0 || !stageAudit.disabled.phase2 ||
          !stageAudit.disabled.phase3 || !stageAudit.disabled.phase4 ||
          !stageAudit.disabled.error) {
        throw new Error(`model ${modelId} stage audit failed: ${JSON.stringify(stageAudit)}`);
      }

      await page.click('[data-mode="phase1"]');
      await page.click('#fitButton');
      await page.waitForTimeout(200);
      const pixels = await page.evaluate(() => window.__primitiveViewerPixelSample());
      if (pixels.uniqueColors < 3 || pixels.nonBackground < 20) {
        throw new Error(`model ${modelId} phase1 is blank: ${JSON.stringify(pixels)}`);
      }
      await page.click('[data-mode="split"]');
      await page.waitForTimeout(100);
      const split = await page.evaluate(() => window.__primitiveViewerDiagnostics());
      if (!split.sourceVisible || split.phase1Meshes < 1) {
        throw new Error(`model ${modelId} split audit failed: ${JSON.stringify(split)}`);
      }
      results.push({ model: modelId, phase1Meshes: diagnostics.phase1Meshes, pixels });
    }
    if (errors.length) throw new Error(`browser errors: ${errors.join(' | ')}`);
    console.log(JSON.stringify({ ok: true, results }));
  } finally {
    await browser.close();
  }
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
