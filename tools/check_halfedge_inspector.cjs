const { chromium } = require('playwright');

const url = process.argv[2];
if (!url) throw new Error('usage: node check_halfedge_inspector.cjs <viewer-url>');

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
    await page.waitForFunction(() => window.__viewerLoadStage === 'complete', null,
      { timeout: 120000 });
    const modes = await page.locator('[data-mode]').evaluateAll(buttons =>
      buttons.map(button => ({
        mode: button.dataset.mode,
        hidden: button.hidden,
        text: button.textContent,
      })));
    for (const hiddenMode of ['phase3', 'phase4', 'error']) {
      if (!modes.find(item => item.mode === hiddenMode)?.hidden)
        throw new Error(`${hiddenMode} should be hidden`);
    }
    await page.fill('#faceInput', '0');
    await page.press('#faceInput', 'Enter');
    const faceStats = await page.locator('#halfedgeStats').innerText();
    if (!faceStats.includes('face_halfedge') || !faceStats.includes('halfedge cycle'))
      throw new Error(`face inspection is incomplete: ${faceStats}`);
    await page.fill('#halfedgeInput', '0');
    await page.click('#inspectHalfedgeButton');
    const halfedgeStats = await page.locator('#halfedgeStats').innerText();
    for (const field of ['origin → destination', 'face', 'next', 'opposite'])
      if (!halfedgeStats.includes(field))
        throw new Error(`halfedge inspection lacks ${field}: ${halfedgeStats}`);
    await page.selectOption('#halfedgeFilter', 'unpaired');
    await page.click('#nextHalfedgeButton');
    const filteredStats = await page.locator('#halfedgeStats').innerText();
    if (!filteredStats.includes('INVALID (unpaired)') &&
        !filteredStats.includes('无匹配 halfedge'))
      throw new Error(`unpaired filter failed: ${filteredStats}`);
    await page.click('[data-mode="phase2"]');
    await page.click('#fitButton');
    await page.screenshot({ path: 'viewer/screenshots/model3-final-90-no-black-material.png',
      fullPage: true });
    const diagnostics = await page.evaluate(() => window.__primitiveViewerDiagnostics());
    if (!diagnostics.phase2Meshes || diagnostics.triangles !== 90)
      throw new Error(`final proxy did not render: ${JSON.stringify(diagnostics)}`);
    const mobile = await browser.newPage({ viewport: { width: 390, height: 844 } });
    mobile.on('pageerror', error => errors.push(error.message));
    await mobile.goto(url, { waitUntil: 'domcontentloaded' });
    await mobile.waitForFunction(() => window.__viewerLoadStage === 'complete', null,
      { timeout: 120000 });
    await mobile.fill('#halfedgeInput', '1');
    await mobile.click('#inspectHalfedgeButton');
    const mobileLayout = await mobile.evaluate(() => {
      const panel = document.querySelector('#panel').getBoundingClientRect();
      return {
        overflow: document.documentElement.scrollWidth > innerWidth,
        panel: { left: panel.left, right: panel.right, top: panel.top, bottom: panel.bottom },
        stats: document.querySelector('#halfedgeStats').innerText,
      };
    });
    if (mobileLayout.overflow || mobileLayout.panel.left < 0 ||
        mobileLayout.panel.right > 390 || mobileLayout.panel.bottom > 844)
      throw new Error(`mobile inspector overflow: ${JSON.stringify(mobileLayout)}`);
    await mobile.screenshot({ path: 'viewer/screenshots/model3-halfedge-inspector-mobile.png',
      fullPage: true });
    await mobile.close();
    if (errors.length) throw new Error(errors.join(' | '));
    console.log(JSON.stringify({ ok: true, modes, faceStats, halfedgeStats,
      filteredStats, diagnostics, mobileLayout }));
  } finally {
    await browser.close();
  }
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
