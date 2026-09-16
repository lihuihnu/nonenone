#!/usr/bin/env python3
"""Download immutable NIST ThermoML inputs and verify their SHA-256 hashes."""

from __future__ import annotations

import hashlib
import json
import pathlib
import urllib.request


SOURCES = {
    "may2015_pr_vle.json": (
        "https://trc.nist.gov/ThermoML/10.1021/acs.jced.5b00610.json",
        "77630e90db70bb6aabfdfa520f61f14cee5076ece0265754140a25f771659662",
    ),
    "messabeb2016_sw_co2_brine.json": (
        "https://trc.nist.gov/ThermoML/10.1021/acs.jced.6b00505.json",
        "71a0762e6ce7c0a89491087cfccee96e507f164b63d7b539ae6d6f60614f715f",
    ),
    "soujanya2010_cpa_meoh_water.json": (
        "https://trc.nist.gov/ThermoML/10.1016/j.jct.2009.11.020.json",
        "8382fec6b97ae7035ef355b3cf9e29beb1cb90e48d809159f8f40b305304d81d",
    ),
    "petropoulou2018_co2_ch4_vle.json": (
        "https://trc.nist.gov/ThermoML/10.1016/j.fluid.2018.01.011.json",
        "f26f16b53fe6c6ad73240215937a01f6e58239c8622a58a9b637d1391d211ac8",
    ),
    "yang2018_co2_ethane_density.json": (
        "https://trc.nist.gov/ThermoML/10.1016/j.jct.2018.02.021.json",
        "2317b9fdfa4bf4892f9f59f3507d2b1f1c48c2e8d21afdae04af5dc008767f30",
    ),
    "frost2014_methane_water_vle.json": (
        "https://trc.nist.gov/ThermoML/10.1021/je400684k.json",
        "a78159d0ab6cfbadebb4028c590dffff5ccdf0acba9e4938bbcbc62379b9fe4a",
    ),
}


def main() -> None:
    destination = pathlib.Path(__file__).resolve().parent / "reference_data"
    destination.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, str]] = []
    for filename, (url, expected_hash) in SOURCES.items():
        path = destination / filename
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "MPMC-SCW-public-benchmark/1.0"},
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read()
        actual_hash = hashlib.sha256(payload).hexdigest()
        if actual_hash != expected_hash:
            raise RuntimeError(
                f"SHA-256 mismatch for {filename}: {actual_hash} != {expected_hash}"
            )
        path.write_bytes(payload)
        manifest.append(
            {
                "file": filename,
                "url": url,
                "sha256": actual_hash,
            }
        )
        print(f"[reference] {filename}: {len(payload)} bytes, sha256={actual_hash}")
    (destination / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
