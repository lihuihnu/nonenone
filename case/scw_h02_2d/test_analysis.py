import importlib.util
import math
import os
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('audit', Path(__file__).with_name('analyze_h02.py'))
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

class AuditTests(unittest.TestCase):
    def test_exact_target_not_later_row(self):
        rows = [dict(pvi=1.98,y=4.),dict(pvi=1.999999997,y=5.),dict(pvi=2.02,y=6.)]
        self.assertLess(abs(audit.interpolate(rows,2.,'y')-5.),2e-7)
    def test_incomplete_coverage_is_not_last_row(self):
        with self.assertRaises(ValueError):
            audit.interpolate([dict(pvi=0.,y=0.),dict(pvi=1.99,y=1.)],2.,'y')
    def test_nan_cannot_satisfy_gate(self):
        with self.assertRaises(ValueError): audit.num(dict(pvi=math.nan),'pvi')
    def test_missing_phase_column_is_not_zero(self):
        with self.assertRaises(KeyError): audit.num({},'cum_Water_OIL_HEAVY_produced_kg')
    def test_no_forward_fill_of_zero_flow_viscosity(self):
        rows=[dict(pvi=0.,mu=math.nan),dict(pvi=0.02,mu=0.002)]
        self.assertTrue(math.isnan(audit.interpolate(rows,0.,'mu')))
        self.assertTrue(math.isnan(audit.interpolate(rows,0.01,'mu')))
    def test_viscosity_reconstruction_endpoints(self):
        for mode in 'BC':
            self.assertAlmostEqual(audit.h02_mu(0.,mode),0.002)
            self.assertAlmostEqual(audit.h02_mu(1.,mode),5e-5)
    def test_nonmonotone_pvi_rejected(self):
        with self.assertRaises(ValueError):
            audit.bracket([dict(pvi=0.),dict(pvi=1.),dict(pvi=0.9)],0.5)
    def test_actual_run30_fails_mass_gate(self):
        location=os.environ.get('H02_RUN30_ROOT')
        if not location: self.skipTest('Set H02_RUN30_ROOT to the unmodified Run-30 full/ artifact')
        with tempfile.TemporaryDirectory() as tmp:
            status=audit.audit(Path(location),Path(tmp))
        self.assertEqual(status['pvi_coverage'],'PASS')
        self.assertEqual(status['numerical_gate'],'FAIL')

if __name__=='__main__': unittest.main()
