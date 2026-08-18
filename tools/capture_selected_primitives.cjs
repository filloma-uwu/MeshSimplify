const path = require('path');
const { chromium } = require('playwright');

const url = process.argv[2];
const outputDirectory = path.resolve(process.argv[3]);
const ids = (process.argv[4] || '1,20').split(',').map(Number);
if (!url || !outputDirectory) throw new Error('usage: ... <url> <output> [ids]');

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
    await page.locator('#otherOpacity').fill('12');
    await page.locator('#otherOpacity').dispatchEvent('input');
    for (const id of ids) {
      await page.selectOption('#primitiveSelect', String(id));
      await page.waitForTimeout(250);
      await page.screenshot({ path: path.join(outputDirectory, `primitive-${id}.png`) });
    }
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
