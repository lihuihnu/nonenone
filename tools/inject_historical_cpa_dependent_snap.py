from pathlib import Path

path = Path('models/include/natural/state/state_codec.hpp')
text = path.read_text()
old = '''            if (std::abs(dependentValue) <=
                NaturalNumerics::dependentCompositionCancellationBoundary)
            {
                dependent += -dependentValue;
            }
'''
new = '''            if (dependentValue < 0.0 &&
                dependentValue >= -NaturalNumerics::phaseEquilibriumTraceComposition)
            {
                dependent += -dependentValue;
            }
'''
if text.count(old) != 1:
    raise SystemExit(f'expected one dependent-snap anchor, found {text.count(old)}')
path.write_text(text.replace(old, new, 1))
print('restored historical dependent-last snapping for CPA regression test')
