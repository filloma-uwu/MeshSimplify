const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');

const url = process.argv[2];
const outputDirectory = path.resolve(process.argv[3] || 'viewer/screenshots/orbit');
if (!url) throw new Error('usage: node capture_orbit_stages.cjs <viewer-url> [output-directory]');

const modes = ['source', 'phase0', 'phase1', 'phase3'];

(async () => {
  fs.mkdirSync(outputDirectory, { recursive: true });
  const browser = await chromium.launch({ headless: true });
  try {
    const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
    await page.goto(url, { waitUntil: 'domcontentloaded' });
    await page.waitForFunction(() => {
      const state = window.__primitiveViewerDiagnostics?.();
      return state?.manifestComplete && !state?.loading;
    }, null, { timeout: 120000 });
    const model = await page.locator('#modelSelect').inputValue();
    const viewport = page.locator('#viewport');
    const box = await viewport.boundingBox();
    const center = { x: box.x + box.width * 0.5, y: box.y + box.height * 0.5 };

    for (let view = 0; view < 8; ++view) {
      for (const mode of modes) {
        const button = page.locator(`[data-mode="${mode}"]`);
        if (await button.isDisabled()) continue;
        await button.click();
        await page.waitForTimeout(100);
        await page.screenshot({
          path: path.join(outputDirectory, `model-${model}-view-${view}-${mode}.png`),
          fullPage: true,
        });
      }
      await page.mouse.move(center.x + 180, center.y);
      await page.mouse.down();
      await page.mouse.move(center.x - 40, center.y, { steps: 12 });
      await page.mouse.up();
      await page.waitForTimeout(120);
    }
  } finally {
    await browser.close();
  }
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
