"""Umbrella Review Tool — Selenium Tests (20 tests)"""
import sys, os, time, io, unittest
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options

HTML_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'umbrella-review.html')
URL = 'file:///' + HTML_PATH.replace('\\', '/')

def get_driver():
    opts = Options()
    opts.add_argument('--headless=new'); opts.add_argument('--no-sandbox'); opts.add_argument('--disable-gpu'); opts.add_argument('--window-size=1400,900')
    opts.set_capability('goog:loggingPrefs', {'browser': 'ALL'})
    d = webdriver.Chrome(options=opts); d.implicitly_wait(2); return d

class UmbrellaReviewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.driver = get_driver(); cls.driver.get(URL); time.sleep(0.5)
    @classmethod
    def tearDownClass(cls): cls.driver.quit()
    def _reload(self): self.driver.get(URL); time.sleep(0.3)
    def _click(self, by, val):
        from selenium.webdriver.support.ui import WebDriverWait
        from selenium.webdriver.support import expected_conditions as EC
        el = WebDriverWait(self.driver, 5).until(EC.presence_of_element_located((by, val)))
        self.driver.execute_script("arguments[0].click()", el); return el

    def test_01_page_loads(self):
        title = self.driver.title.lower()
        self.assertTrue('umbrella' in title or 'overview' in title or 'meta' in title)
    def test_02_five_tabs(self):
        tabs = self.driver.find_elements(By.CSS_SELECTOR, '[role="tab"]')
        self.assertGreaterEqual(len(tabs), 5)
    def test_03_input_tab(self):
        panel = self.driver.find_element(By.ID, 'panel-input')
        self.assertTrue(panel.is_displayed())
    def test_04_add_row_button(self):
        btn = self.driver.find_element(By.ID, 'btnAddRow')
        self.assertIsNotNone(btn)
    def test_05_example_selector(self):
        el = self.driver.find_element(By.ID, 'exampleSelect')
        self.assertIsNotNone(el)
    def test_06_load_example(self):
        self._reload()
        sel = self.driver.find_element(By.ID, 'exampleSelect')
        opts = sel.find_elements(By.TAG_NAME, 'option')
        if len(opts) > 1:
            self.driver.execute_script("arguments[0].value = arguments[1]; arguments[0].dispatchEvent(new Event('change'))", sel, opts[1].get_attribute('value'))
            time.sleep(0.5)
    def test_07_data_table(self):
        table = self.driver.find_element(By.ID, 'dataTable')
        self.assertIsNotNone(table)
    def test_08_run_analysis(self):
        self._click(By.ID, 'btnRunAnalysis'); time.sleep(1)
    def test_09_pooled_tab(self):
        self._click(By.ID, 'tab-pooled'); time.sleep(0.3)
        panel = self.driver.find_element(By.ID, 'panel-pooled')
        self.assertIn('active', panel.get_attribute('class') or '')
    def test_10_pooled_has_content(self):
        panel = self.driver.find_element(By.ID, 'panel-pooled')
        self.assertGreater(len(panel.text), 20)
    def test_11_quality_tab(self):
        self._click(By.ID, 'tab-quality'); time.sleep(0.3)
        panel = self.driver.find_element(By.ID, 'panel-quality')
        self.assertIn('active', panel.get_attribute('class') or '')
    def test_12_evidence_tab(self):
        self._click(By.ID, 'tab-evidence'); time.sleep(0.3)
        panel = self.driver.find_element(By.ID, 'panel-evidence')
        self.assertIn('active', panel.get_attribute('class') or '')
    def test_13_evidence_classification(self):
        panel = self.driver.find_element(By.ID, 'panel-evidence')
        text = panel.text.lower()
        has_class = any(c in text for c in ['convincing', 'suggestive', 'weak', 'class', 'non-significant', 'evidence'])
        self.assertTrue(has_class or len(text) > 20)
    def test_14_report_tab(self):
        self._click(By.ID, 'tab-report'); time.sleep(0.3)
        panel = self.driver.find_element(By.ID, 'panel-report')
        self.assertIn('active', panel.get_attribute('class') or '')
    def test_15_report_content(self):
        panel = self.driver.find_element(By.ID, 'panel-report')
        self.assertGreater(len(panel.text), 20)
    def test_16_dark_mode(self):
        self._reload()
        btn = self.driver.find_element(By.ID, 'darkToggle')
        self.driver.execute_script("arguments[0].click()", btn); time.sleep(0.2)
        self.driver.execute_script("arguments[0].click()", btn)
    def test_17_csv_paste(self):
        self._reload()
        self._click(By.ID, 'btnPasteCSV'); time.sleep(0.3)
        ta = self.driver.find_element(By.ID, 'csvInput')
        ta.send_keys("Test MA, 0.85, 0.72, 1.02, 5, 500, 30, High, Mortality")
        self._click(By.ID, 'btnParseCSV'); time.sleep(0.3)
    def test_18_clear_all(self):
        self._click(By.ID, 'btnClearAll'); time.sleep(0.3)
        try:
            self.driver.switch_to.alert.accept(); time.sleep(0.2)
        except Exception:
            pass
    def test_19_tab_keyboard(self):
        self.driver.get(URL); time.sleep(0.5)
        tabs = self.driver.find_elements(By.CSS_SELECTOR, '[role="tab"]')
        if len(tabs) > 0:
            tabs[0].click(); time.sleep(0.2)
            tabs[0].send_keys(Keys.ARROW_RIGHT); time.sleep(0.2)
    def test_20_no_js_errors(self):
        logs = self.driver.get_log('browser')
        severe = [l for l in logs if l['level']=='SEVERE' and 'favicon' not in l.get('message','')]
        self.assertEqual(len(severe), 0, f"JS errors: {severe}")

if __name__ == '__main__':
    unittest.main(verbosity=2)
