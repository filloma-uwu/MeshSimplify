const { test, expect } = require('playwright/test');

const viewerUrl = process.env.PRIMITIVE_VIEWER_URL ||
  'http://127.0.0.1:8091/viewer/primitive_analysis.html' +
  '?manifest=/outputs/spatial_group_gap_validation/viewer_manifest.json&model=3';

function expectVectorUnchanged(current, initial, tolerance = 1e-8) {
  expect(current).toHaveLength(initial.length);
  current.forEach((value, index) => expect(Math.abs(value - initial[index])).toBeLessThan(tolerance));
}

test('proxy OBJ modes preserve the camera and expose no legacy responsibility UI', async ({ page }) => {
  test.setTimeout(120_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto(viewerUrl);
  await page.waitForFunction(() => {
    const state = window.__primitiveViewerDiagnostics?.();
    return state?.model === 3 && state?.manifestComplete && !state?.loading;
  });

  const initial = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  expect(initial.primitiveCount).toBeGreaterThan(0);
  expect(initial.primitiveTypeCounts.frustum || 0).toBe(0);
  expect(initial.primitiveTypeCounts.obb).toBeUndefined();

  for (const mode of ['source', 'phase1', 'phase2', 'phase3', 'phase4', 'split']) {
    await page.locator(`[data-mode="${mode}"]`).click();
    await page.waitForFunction(expected => window.__primitiveViewerDiagnostics?.().mode === expected, mode);
    const current = await page.evaluate(() => window.__primitiveViewerDiagnostics());
    expectVectorUnchanged(current.cameraPosition, initial.cameraPosition);
    expectVectorUnchanged(current.cameraTarget, initial.cameraTarget);
  }

  await expect(page.locator('[data-mode="regions"]')).toHaveCount(0);
  await expect(page.locator('#strength')).toHaveCount(0);
  await expect(page.locator('#stats')).toContainText('图元数');
  await expect(page.locator('#stats')).not.toContainText('责任');
  await expect(page.locator('input[data-type="disk"]')).toHaveCount(1);
  await expect(page.locator('input[data-type="annulus"]')).toHaveCount(1);
  await expect(page.locator('input[data-type="cylindricalband"]')).toHaveCount(1);
  await expect(page.locator('input[data-type="conicalband"]')).toHaveCount(1);
});

test('batched models keep primitive selection and reset selection opacity on model change', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto(
    'http://127.0.0.1:8093/viewer/primitive_analysis.html' +
    '?manifest=/outputs/surface_primitives_v107_official_scope/viewer_manifest.json&model=12',
  );
  await page.waitForFunction(() => {
    const state = window.__primitiveViewerDiagnostics?.();
    return state?.model === 12 && !state?.loading;
  });

  let state = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  expect(state.primitiveOptionCount).toBe(state.primitiveCount);
  expect(state.primitiveOptionCount).toBeGreaterThan(5000);

  await page.locator('#primitiveSelect').selectOption('8072');
  await page.waitForFunction(() =>
    window.__primitiveViewerDiagnostics?.().selectedPrimitive === '8072');
  state = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  expect(state.selectionOverlayCount).toBe(2);
  expect(state.primitiveOpacity).toBeCloseTo(0.15, 5);
  expect(state.triangulatedOpacity).toBeCloseTo(0.15, 5);

  await page.route('**/models/13/model.json', async route => {
    await new Promise(resolve => setTimeout(resolve, 750));
    await route.continue();
  });

  await page.locator('#modelSelect').selectOption('13');
  await page.waitForFunction(() => {
    const current = window.__primitiveViewerDiagnostics?.();
    return current?.selectedPrimitive === 'all' && current?.loading;
  });
  state = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  expect(state.selectionOverlayCount).toBe(0);
  expect(state.allBaseMaterialsOpaque).toBe(true);

  await page.waitForFunction(() => {
    const current = window.__primitiveViewerDiagnostics?.();
    return current?.model === 13 && !current?.loading;
  });
  state = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  expect(state.selectedPrimitive).toBe('all');
  expect(state.primitiveOptionCount).toBe(state.primitiveCount);
  expect(state.primitiveOpacity).toBeCloseTo(1, 5);
  expect(state.triangulatedOpacity).toBeCloseTo(1, 5);
  expect(state.allBaseMaterialsOpaque).toBe(true);
  expect(state.selectionOverlayCount).toBe(0);
});

test('every current batch model exposes selectable primitives and restores opacity', async ({ page }) => {
  test.setTimeout(300_000);
  const modelIds = [2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20];
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto(
    'http://127.0.0.1:8093/viewer/primitive_analysis.html' +
    '?manifest=/outputs/spatial_group_gap_validation/viewer_manifest.json&model=2',
  );

  for (const modelId of modelIds) {
    if (modelId !== modelIds[0]) {
      await page.locator('#modelSelect').selectOption(String(modelId));
    }
    await page.waitForFunction(expectedModel => {
      const current = window.__primitiveViewerDiagnostics?.();
      return current?.model === expectedModel && !current?.loading;
    }, modelId);

    let state = await page.evaluate(() => window.__primitiveViewerDiagnostics());
    expect(state.primitiveCount).toBeGreaterThan(0);
    expect(state.primitiveOptionCount).toBe(state.primitiveCount);
    expect(state.selectedPrimitive).toBe('all');
    expect(state.allBaseMaterialsOpaque).toBe(true);

    const lastPrimitive = String(state.primitiveCount - 1);
    await page.locator('#primitiveSelect').selectOption(lastPrimitive);
    await page.waitForFunction(expectedPrimitive =>
      window.__primitiveViewerDiagnostics?.().selectedPrimitive === expectedPrimitive,
    lastPrimitive);
    state = await page.evaluate(() => window.__primitiveViewerDiagnostics());
    if (state.selectionOverlayCount > 0) {
      expect(state.selectionOverlayCount).toBe(2);
      expect(state.primitiveOpacity).toBeCloseTo(0.15, 5);
      expect(state.triangulatedOpacity).toBeCloseTo(0.15, 5);
    } else {
      expect(state.selectedBaseMeshCount).toBe(2);
      expect(state.selectedBaseMaterialsOpaque).toBe(true);
      expect(state.nonSelectedBaseMaterialsFaded).toBe(true);
    }

    await page.locator('#primitiveSelect').selectOption('all');
    await page.waitForFunction(() => {
      const current = window.__primitiveViewerDiagnostics?.();
      return current?.selectedPrimitive === 'all' &&
        current?.allBaseMaterialsOpaque === true;
    });
  }
});

test('open-error mode preserves the camera and locates the maximum-distance pair', async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto(
    'http://127.0.0.1:8093/viewer/primitive_analysis.html' +
    '?manifest=/outputs/spatial_group_gap_validation/viewer_manifest.json&model=2',
  );
  await page.waitForFunction(() => {
    const state = window.__primitiveViewerDiagnostics?.();
    return state?.model === 2 && !state?.loading;
  });

  const initial = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  await page.locator('[data-mode="error"]').click();
  const error = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  expectVectorUnchanged(error.cameraPosition, initial.cameraPosition);
  expectVectorUnchanged(error.cameraTarget, initial.cameraTarget);
  expect(error.errorVisible).toBe(true);
  expect(error.errorBoundaryPoints).toBeGreaterThanOrEqual(0);
  expect(error.maximumErrorDistance).toBeGreaterThan(0);
  expect(error.maximumErrorSegmentCount).toBe(1);
  await expect(page.locator('#stats')).toContainText('严格包围验证');
  await expect(page.locator('#stats')).toContainText('最大误差');

  await page.locator('#focusErrorButton').click();
  const focused = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  expect(focused.cameraTarget).not.toEqual(initial.cameraTarget);
});

test('maximum-error regeneration reloads the proxy without resetting the camera', async ({ page }) => {
  test.setTimeout(180_000);
  let regenerationRequest = null;
  page.on('request', request => {
    if (request.url().endsWith('/api/regenerate'))
      regenerationRequest = request.postDataJSON();
  });
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto(
    'http://127.0.0.1:8094/viewer/primitive_analysis.html' +
    '?manifest=/outputs/model3_max_error_workload_only/viewer_manifest.json&model=3',
  );
  await page.waitForFunction(() => {
    const state = window.__primitiveViewerDiagnostics?.();
    return state?.model === 3 && !state?.loading;
  });

  const initial = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  await page.locator('#maximumErrorInput').fill('80');
  await page.locator('#regenerateButton').click();
  await page.waitForFunction(() => {
    const state = window.__primitiveViewerDiagnostics?.();
    return !state?.loading && state?.maximumErrorLimit === 80 &&
      state?.regenerateStatus.includes('80');
  }, null, { timeout: 150_000 });

  const regenerated = await page.evaluate(() => window.__primitiveViewerDiagnostics());
  expectVectorUnchanged(regenerated.cameraPosition, initial.cameraPosition);
  expectVectorUnchanged(regenerated.cameraTarget, initial.cameraTarget);
  expect(regenerated.maximumErrorDistance).toBeLessThanOrEqual(80 + 1e-8);
  expect(regenerated.maximumErrorInput).toBe(80);
  expect(regenerationRequest.manifest).toBe(
    '/outputs/model3_max_error_workload_only/viewer_manifest.json');
});
