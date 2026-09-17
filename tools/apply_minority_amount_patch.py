#!/usr/bin/env python3
from pathlib import Path
import base64, zlib, subprocess

PATCH_B64 = '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'
patch = zlib.decompress(base64.b64decode(PATCH_B64))
Path('/tmp/minority_amount.patch').write_bytes(patch)
subprocess.run(['git','apply','--check','/tmp/minority_amount.patch'], check=True)
subprocess.run(['git','apply','/tmp/minority_amount.patch'], check=True)

# Patch the live current SW regression instead of replacing it with the older handoff copy.
from pathlib import Path

path = Path('test/src/unit/sw_phase_ordering_test.cpp')
s = path.read_text()

def replace_once(old: str, new: str, label: str) -> None:
    global s
    count = s.count(old)
    if count != 1:
        raise RuntimeError(f'{label}: expected exactly one match, found {count}')
    s = s.replace(old, new, 1)

replace_once(
'''void writeComposition(
    std::array<double, Indices::numPrimaryVariables>& primary,
    const std::array<int, Indices::numIndependentCompositionsPerPhase>& indices,
    const Composition& composition)
{
    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        primary[static_cast<std::size_t>(
            indices[static_cast<std::size_t>(component)])] =
            composition[static_cast<std::size_t>(component)];
    }
}
''',
'''void writeComposition(
    std::array<double, Indices::numPrimaryVariables>& primary,
    MPMC::CompositionalPhase phase,
    const Composition& composition)
{
    const auto &indices = phase == MPMC::CompositionalPhase::Oil
        ? Indices::Primary::liquidComposition
        : (phase == MPMC::CompositionalPhase::Gas
               ? Indices::Primary::vaporComposition
               : Indices::Primary::waterComposition);
    const double amountScale = phase == MPMC::CompositionalPhase::Oil
        ? primary[Indices::Primary::liquidSaturation]
        : 1.0;
    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        const auto c = static_cast<std::size_t>(component);
        primary[static_cast<std::size_t>(indices[c])] =
            amountScale * composition[c];
    }
}

Composition readComposition(
    const std::array<double, Indices::numPrimaryVariables>& primary,
    MPMC::CompositionalPhase phase)
{
    const auto &indices = phase == MPMC::CompositionalPhase::Oil
        ? Indices::Primary::liquidComposition
        : (phase == MPMC::CompositionalPhase::Gas
               ? Indices::Primary::vaporComposition
               : Indices::Primary::waterComposition);
    const double amountScale = phase == MPMC::CompositionalPhase::Oil
        ? primary[Indices::Primary::liquidSaturation]
        : 1.0;
    if (!(amountScale > 0.0))
        throw std::runtime_error("active Oil composition requires positive saturation");

    Composition composition{};
    composition.back() = 1.0;
    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        const auto c = static_cast<std::size_t>(component);
        composition[c] = primary[static_cast<std::size_t>(indices[c])] / amountScale;
        composition.back() -= composition[c];
    }
    return composition;
}
''',
'composition helper')

replace_once(
'''    writeComposition(primary, Indices::Primary::liquidComposition, co2Rich);
    writeComposition(primary, Indices::Primary::vaporComposition, nc10Rich);
    writeComposition(primary, Indices::Primary::waterComposition, waterRich);
''',
'''    writeComposition(primary, MPMC::CompositionalPhase::Oil, co2Rich);
    writeComposition(primary, MPMC::CompositionalPhase::Gas, nc10Rich);
    writeComposition(primary, MPMC::CompositionalPhase::Water, waterRich);
''',
'canonicalization writes')

replace_once(
'''    require(
        primary[Indices::Primary::vaporComposition[1]] >
            primary[Indices::Primary::liquidComposition[1]],
        "SW Newton state must map the CO2-rich phase to gas");
    require(
        primary[Indices::Primary::liquidComposition[0]] <
            primary[Indices::Primary::vaporComposition[0]],
        "SW Newton state must map the lower-water nonaqueous phase to oil");
''',
'''    const auto canonicalOil =
        readComposition(primary, MPMC::CompositionalPhase::Oil);
    const auto canonicalGas =
        readComposition(primary, MPMC::CompositionalPhase::Gas);
    require(canonicalGas[1] > canonicalOil[1],
            "SW Newton state must map the CO2-rich phase to gas");
    require(canonicalOil[0] < canonicalGas[0],
            "SW Newton state must map the lower-water nonaqueous phase to oil");
''',
'canonicalization assertions')

replace_once(
'''    writeComposition(primary, Indices::Primary::liquidComposition, oilIterate);
    writeComposition(primary, Indices::Primary::vaporComposition, overall);
    writeComposition(primary, Indices::Primary::waterComposition, waterIterate);
''',
'''    writeComposition(primary, MPMC::CompositionalPhase::Oil, oilIterate);
    writeComposition(primary, MPMC::CompositionalPhase::Gas, overall);
    writeComposition(primary, MPMC::CompositionalPhase::Water, waterIterate);
''',
'transient writes')

replace_once(
'''        require(std::abs(
                    primary[static_cast<std::size_t>(
                        Indices::Primary::liquidComposition[c])] -
                    oilIterate[c]) < 1.0e-14,
                "SW stability check must not overwrite oil Newton iterate");
''',
'''        const auto decodedOil =
            readComposition(primary, MPMC::CompositionalPhase::Oil);
        require(std::abs(decodedOil[c] - oilIterate[c]) < 1.0e-14,
                "SW stability check must not overwrite oil Newton iterate");
''',
'transient oil assertion')

path.write_text(s)
print('patched current SW phase-ordering regression for qO coordinates')
