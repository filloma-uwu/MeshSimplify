const path = require('path');
const { chromium } = require('playwright');

const url = process.argv[2];
const outputDirectory = path.resolve(process.argv[3]);
if (!url || !outputDirectory) {
  throw new Error('usage: node capture_primitive_front_back.cjs <url> <output-directory>');
}

(async () => {
  const browser = await chromium.launch({ headless: true });
  try {
    const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
    await page.goto(url, { waitUntil: 'domcontentloaded' });
    await page.waitForFunction(() => {
      const state = window.__primitiveViewerDiagnostics?.();
      return state?.manifestComplete && !state?.loading;
    }, null, { timeout: 120000 });
    await page.click('[data-mode="phase4"]');
    await page.click('#fitButton');
    await page.waitForTimeout(300);
    await page.screenshot({ path: path.join(outputDirectory, 'front.png') });

    const canvas = page.locator('canvas');
    const bounds = await canvas.boundingBox();
    if (!bounds) throw new Error('viewer canvas is not visible');
    const x = bounds.x + bounds.width * 0.5;
    const y = bounds.y + bounds.height * 0.5;
    await page.mouse.move(x, y);
    await page.mouse.down({ button: 'left' });
    // OrbitControls maps a horizontal drag of half the canvas height to pi.
    await page.mouse.move(x + bounds.height * 0.5, y, { steps: 40 });
    await page.mouse.up({ button: 'left' });
    await page.waitForTimeout(300);
    await page.screenshot({ path: path.join(outputDirectory, 'back.png') });
  } finally {
    await browser.close();
  }
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
