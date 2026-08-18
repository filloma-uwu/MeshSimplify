const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');

const url = process.argv[2];
const outputDirectory = path.resolve(process.argv[3] || 'viewer/screenshots/stages');
const [viewportWidth, viewportHeight] = (process.env.PQSS_VIEWPORT || '1440x900')
  .split('x').map(Number);
if (!url) {
  throw new Error('usage: node capture_primitive_stages.cjs <viewer-url> [output-directory]');
}
if (!(viewportWidth > 0) || !(viewportHeight > 0)) {
  throw new Error('PQSS_VIEWPORT must use WIDTHxHEIGHT, for example 1158x919');
}

const modes = ['source', 'phase1', 'phase2', 'phase3', 'phase4', 'error', 'split'];

(async () => {
  fs.mkdirSync(outputDirectory, { recursive: true });
  const browser = await chromium.launch({ headless: true });
  try {
    const page = await browser.newPage({
      viewport: { width: viewportWidth, height: viewportHeight },
    });
    await page.goto(url, { waitUntil: 'domcontentloaded' });
    await page.waitForFunction(() => {
      const state = window.__primitiveViewerDiagnostics?.();
      return state?.manifestComplete && !state?.loading;
    }, null, { timeout: 120000 });
    const model = await page.locator('#modelSelect').inputValue();
    for (const mode of modes) {
      const button = page.locator(`[data-mode="${mode}"]`);
      if (await button.isDisabled()) continue;
      await button.click();
      await page.click('#fitButton');
      await page.waitForTimeout(150);
      await page.screenshot({
        path: path.join(outputDirectory, `model-${model}-${mode}.png`),
        fullPage: true,
      });
    }
  } finally {
    await browser.close();
  }
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
