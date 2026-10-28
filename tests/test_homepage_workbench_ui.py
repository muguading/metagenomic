from __future__ import annotations

import shutil
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class HomePageWorkbenchTests(unittest.TestCase):
    def test_homepage_is_list_first_and_keeps_workflow_hooks(self):
        source = (ROOT / "bac_analysis_portal/templates/index.html").read_text()
        self.assertTrue("filename='css/domains/submission-home.css'" in source, "Required homepage markup missing")
        self.assertTrue('class="workstation-ui"' in source, "Required homepage markup missing")
        self.assertTrue('id="submission-module-guide"' in source, "Required homepage markup missing")
        self.assertTrue('role="radiogroup"' in source, "Required homepage markup missing")
        self.assertTrue('新建分析任务' in source, "Required homepage markup missing")
        self.assertFalse('Task Intake' in source, "Redundant homepage markup remains")
        self.assertFalse('Submission Summary' in source, "Redundant homepage markup remains")
        self.assertFalse('这里以后可以继续扩展更多模块' in source, "Redundant homepage markup remains")
        for key in ('bacteria', 'virus', 'metagenome', 'community', 'pathosource'):
            self.assertTrue(f'id="module-choice-{key}"' in source, "Required homepage markup missing")
            self.assertRegex(source, rf'type="radio" name="analysis_module" value="{key}"')
        for step in range(1, 5):
            self.assertTrue(f'id="submission-step-{step}"' in source, "Required homepage markup missing")
            self.assertTrue(f'data-submission-band-step="{step}"' in source, "Required homepage markup missing")
        for target in ('submission-tab', 'queue-tab', 'server-tab', 'host-tab', 'database-tab',
                       'modeling-tab', 'project-tab', 'audit-tab', 'admin-tab', 'database-extra-tab'):
            self.assertEqual(source.count(f'data-tab-target="{target}"'), 1)

    def test_new_styles_keep_focus_and_narrow_screen_access(self):
        source = (ROOT / 'bac_analysis_portal/static/css/domains/submission-home.css').read_text()
        for selector in (':focus-visible', ':focus-within', '@media', '.sidebar-collapsed',
                         '.module-choice-card.hidden', 'minmax(0, 1fr)', 'prefers-reduced-motion'):
            self.assertTrue(selector in source, "Required homepage markup missing")
        self.assertFalse('linear-gradient' in source, "Redundant homepage markup remains")
        self.assertFalse('backdrop-filter: blur(' in source, "Redundant homepage markup remains")

    def test_step_runtime_updates_context_and_current_step(self):
        source = (ROOT / 'bac_analysis_portal/static/app.js').read_text()
        self.assertTrue('renderSubmissionModuleGuide' in source, "Required homepage markup missing")
        self.assertTrue('setAttribute("aria-current", "step")' in source, "Required homepage markup missing")
        self.assertTrue('setAttribute("data-submission-step", String(currentStep))' in source, "Required homepage markup missing")
        self.assertTrue('classList.toggle("hidden", currentStep === 1)' in source, "Required homepage markup missing")

    @unittest.skipUnless(shutil.which('node'), 'Node.js required for frontend module execution')
    def test_module_guide_renders_each_module_and_no_false_completion(self):
        path = ROOT / 'bac_analysis_portal/static/workbench/domains/submission/home.js'
        script = r"""
          const fs = require('fs'); const assert = require('node:assert/strict');
          (async () => {
            const source = fs.readFileSync(process.argv[1], 'utf8');
            const mod = await import('data:text/javascript;base64,'+Buffer.from(source).toString('base64'));
            for(const key of ['bacteria','virus','metagenome','community','pathosource']) {
              const nodes = Object.fromEntries(['title','description','steps','outputs'].map(k => [k, {textContent:'',innerHTML:''}]));
              const root = {hidden:false, querySelector: selector => nodes[selector.match(/data-guide-(.*)\]/)[1]]};
              const guide = mod.renderSubmissionModuleGuide(root, key);
              assert.equal(root.hidden, false); assert.ok(guide.label.endsWith('分析'));
              assert.equal(nodes.title.textContent, guide.label + ' · 流程说明');
              assert.ok(nodes.steps.innerHTML.includes('<li>'));
              assert.ok(nodes.outputs.innerHTML.includes('<li>'));
              assert.ok(!nodes.steps.innerHTML.includes('is-complete'));
              assert.equal(mod.getSubmissionModuleGuide(key), guide);
            }
            const missing = {hidden:false};
            assert.equal(mod.renderSubmissionModuleGuide(missing, 'not-a-module'), null);
            assert.equal(missing.hidden, true);
            console.log('PASS: five module guides, unknown module, no false completion');
          })().catch(error => {console.error(error); process.exit(1);});
        """
        result = subprocess.run(['node', '-e', script, str(path)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()
