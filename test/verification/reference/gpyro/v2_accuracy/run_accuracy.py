#!/usr/bin/env python3
"""Run/extract the seven pinned Gpyro V2 accuracy cases; never collect timings.

Usage: python3 run_accuracy.py GPYRO WORKSPACE [--extract-only] [--compiler TEXT]
WORKSPACE contains runs/<case> and refined/<case>. Results are written beside
this script. Each refined run halves dt and (for transport) the grid spacing.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import math
import os
from pathlib import Path
import re
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
REVISION = 'bff0126c68f2a35c084ef04490ae327185f281a3'
CASES = ('vgt07a', 'vgt07b', 'vgt08a', 'vgt08b', 'vrx05a', 'vrx05b', 'vrx05c')
R = 8.31446261815324
SUCCESS = 'End of simulation reached successfully. Shutting down.'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, names, rows):
    with path.open('w', newline='') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(names)
        writer.writerows(rows)


def numeric_rows(path, skip=0):
    with path.open() as f:
        rows = list(csv.reader(f))[skip:]
    result = [[float(x) for x in row if x.strip()] for row in rows]
    assert result and len({len(row) for row in result}) == 1, path
    assert all(math.isfinite(v) for row in result for v in row), path
    return result


def erfcx(x):
    if x < 26:
        return math.exp(x*x)*math.erfc(x)
    term = total = 1.0
    for n in range(1, 20):
        term *= -(2*n-1)/(2*x*x)
        total += term
    return total/(math.sqrt(math.pi)*x)


def diffusion_exact(x, t, hm):
    if t == 0:
        return 0.0
    rho = 1e5*.032/(R*300)
    b = hm/(.6*rho*2e-5)
    u = x/(2*math.sqrt(2e-5*t))
    return .22*(math.erfc(u)-math.exp(-u*u)*erfcx(u+b*math.sqrt(2e-5*t)))


def reaction_exact(t, heat, sy):
    a = math.exp(-.2*t)
    cp = sy*1900+(1-sy)*1000
    temp = 1000+heat/(1200-cp)*math.log1p((cp-1200)*(1-a)/1200)
    masses = (1000*a, 1000*sy*(1-a), 1000*(1-sy)*(1-a))
    enthalpy = .001*sum(c*m for c, m in zip((1200, 1900, 1000), masses))*(temp-298.15)
    return temp, masses, enthalpy


def extract(path, case, destination):
    summary = numeric_rows(path/f'{case}_summary_01_0001.csv', skip=1)
    target_end = 400 if case.startswith('vgt07') else 1 if case.startswith('vgt08') else 20
    assert abs(summary[-1][0]-target_end) < 1e-8, (case, summary[-1][0])
    metrics = {}
    if case.startswith('vgt07'):
        hm = 1e-4 if case.endswith('a') else 1e-2
        profiles = numeric_rows(path/f'{case}_profile_01_01_YJ(01).csv')
        x = profiles[0][1:]
        assert abs(x[0]) < 1e-12 and abs(x[-1]-1) < 1e-12
        rho = 1e5*32/(8314*300)  # RHOGOFT's native constant, not the benchmark R.
        b = hm/(.6*(1e5*.032/(R*300))*2e-5)
        rows, history = [], []
        selected = [row for row in profiles[1:] if row[0] in (0., 50., 150., 400.)]
        assert len(selected) == 4
        profile_error = flux_error = inventory_error = 0.
        for t, *y in selected:
            assert len(y) == len(x) and all(-1e-9 <= v <= .22+1e-9 for v in y)
            assert abs(y[-1])*2**-23 < 3e-22, 'non-negligible back-film leakage'
            flux = hm*(.22-y[0])
            mass = .6*rho*sum((x[i+1]-x[i])*(y[i+1]+y[i])/2 for i in range(len(x)-1))
            history.append((t, flux, mass))
            if t > 0:
                rows.extend((t, xx, yy) for xx, yy in zip(x, y))
                profile_error = max(profile_error, max(abs(yy-diffusion_exact(xx, t, hm))/.22 for xx, yy in zip(x, y)))
                flux_exact = hm*.22*erfcx(b*math.sqrt(2e-5*t))
                mass_exact = .6*(1e5*.032/(R*300))*.22*(2*math.sqrt(2e-5*t/math.pi)+(erfcx(b*math.sqrt(2e-5*t))-1)/b)
                flux_error = max(flux_error, abs(flux/flux_exact-1))
                inventory_error = max(inventory_error, abs(mass/mass_exact-1))
        drows = numeric_rows(path/f'{case}_profile_01_02_D12.csv')[2:]
        assert all(abs(d/2e-5-1) < 1e-8 for row in drows for d in row[1:])
        assert max(abs(row[1]+273.15-300) for row in summary) < 1e-6
        write_csv(destination/f'{case}_profiles.csv', ['t_s', 'depth_m', 'Y'], rows)
        write_csv(destination/f'{case}_history.csv', ['t_s', 'influx_kg_m2_s', 'inventory_kg_m2'], history)
        metrics = dict(profile_over_Yinf=profile_error, relative_influx=flux_error, relative_inventory=inventory_error)
    elif case.startswith('vgt08'):
        profiles = numeric_rows(path/f'{case}_profile_01_01_PRESSURE.csv')
        x = profiles[0][1:]
        assert abs(x[0]) < 1e-12 and abs(x[-1]-.1) < 1e-12
        row = next(row for row in profiles[1:] if row[0] == 1.)
        slope = 2*2e-5*R*300*10/(5e-10*.029)
        viscosity = numeric_rows(path/f'{case}_profile_01_03_D12.csv')
        for pressure_row, viscosity_row in zip(profiles[2:], viscosity[2:]):
            assert pressure_row[0] == viscosity_row[0]
            assert all(abs(d*(p+1e5)/(2e-5*R*300/.029)-1) < 1e-8
                       for p, d in zip(pressure_row[1:], viscosity_row[1:]))
        rows = [(1., .1-xx, pp+1e5) for xx, pp in zip(x, row[1:])]
        exact = lambda z: math.sqrt(1e10+slope*(.1-z if case.endswith('a') else z))
        history = []
        for t, temp_c, net_g, top_z_g, *_ in summary:
            if t == 0:
                continue  # The imposed inlet is not initialized in Gpyro's t=0 dump.
            top = -top_z_g/1000
            bottom = net_g/1000-top  # MLR integrates both outward boundary fluxes.
            inlet, outlet = (bottom, top) if case.endswith('a') else (top, bottom)
            history.append((t, inlet, outlet, inlet+outlet))
        assert max(abs(row[1]+273.15-300) for row in summary) < 1e-6
        write_csv(destination/f'{case}_profiles.csv', ['t_s', 'z_m', 'P_Pa'], sorted(rows))
        write_csv(destination/f'{case}_history.csv', ['t_s', 'inlet_kg_m2_s', 'outlet_kg_m2_s', 'net_kg_m2_s'], history)
        metrics = dict(P2_over_rise=max(abs(p*p-exact(z)**2)/(slope*.1) for _, z, p in rows),
                       pressure_Pa=max(abs(p-exact(z)) for _, z, p in rows),
                       relative_flux=max(max(abs(i+10), abs(o-10), abs(net))/10 for _, i, o, net in history))
    else:
        heat = {'a': 1e6, 'b': -1e5, 'c': 0}[case[-1]]
        sy = .4 if case.endswith('c') else 1.
        rows = []
        metrics = dict(temperature_K=0., species_over_initial=0., enthalpy_over_initial=0., total_mass_over_initial=0.)
        for t, tc, a, b, gas, *_ in summary:
            temp = tc+273.15
            enthalpy = .001*(1200*a+1900*b+1000*gas)*(temp-298.15)
            te, me, he = reaction_exact(t, heat, sy)
            metrics['temperature_K'] = max(metrics['temperature_K'], abs(temp-te))
            metrics['species_over_initial'] = max(metrics['species_over_initial'], *(abs(v-e)/1000 for v, e in zip((a,b,gas), me)))
            metrics['enthalpy_over_initial'] = max(metrics['enthalpy_over_initial'], abs(enthalpy-he)/(1200*(1000-298.15)))
            metrics['total_mass_over_initial'] = max(metrics['total_mass_over_initial'], abs((a+b+gas)/1000-1))
            rows.append((t,temp,a,b,gas,enthalpy))
        write_csv(destination/f'{case}_history.csv', ['t_s','T_K','A_kg_m3','B_kg_m3','G_kg_m3','H_J_m2'], rows)
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('workspace', type=Path)
    parser.add_argument('--extract-only', action='store_true')
    parser.add_argument('--compiler', default='unspecified')
    args = parser.parse_args()
    exe = args.executable.resolve()
    work = args.workspace.resolve()
    import json
    provenance = dict(source_repository='https://github.com/Gedeben/GPYRO-V2', source_revision=REVISION,
                      executable_sha256=digest(exe), compiler=args.compiler,
                      generated_utc=datetime.now(timezone.utc).isoformat(), OMP_NUM_THREADS=1,
                      purpose='accuracy only; no timing measurements retained', runs=[])
    metrics_rows = []
    # Time-only triplet at the base mesh, then space-only refinement at dt/4.
    # Keep the original coupled refinement as an explicitly separate sample.
    levels = {'runs': (1, 1), 'refined': (2, 2), 'time_half': (1, 2),
              'time_quarter': (1, 4), 'space_half': (2, 4)}
    for level, (space_factor, time_factor) in levels.items():
        for case in CASES:
            text = (HERE/'decks'/f'{case}.data').read_text()
            if level != 'runs':
                dt = float(re.search(r'DT0=([^,]+)', text)[1])
                text = re.sub(r'DT0=[^,]+', f'DT0={dt/time_factor:.16g}', text)
                if case.startswith('vgt'):
                    n = int(re.search(r'NCELLZ\(1\)=(\d+)', text)[1])
                    text = re.sub(r'NCELLZ\(1\)=\d+', f'NCELLZ(1)={space_factor*(n-1)+1}', text)
            path = work/level/case
            path.mkdir(parents=True, exist_ok=True)
            deck = path/f'{case}.data'
            if args.extract_only:
                assert deck.read_text() == text, f'deck mismatch: {deck}'
            else:
                deck.write_text(text)
                with (path/'run.log').open('w') as log:
                    subprocess.run([str(exe), deck.name], cwd=path, env={**os.environ, 'OMP_NUM_THREADS':'1'},
                                   stdout=log, stderr=subprocess.STDOUT, check=True, timeout=240)
            log = (path/'run.log').read_text()
            assert SUCCESS in log and REVISION in log, f'incomplete or wrong-version run: {path}'
            assert not re.search(r'NaN|Infinity|Failed timestep', log), f'nonfinite/rejected step: {path}'
            dest = HERE if level == 'runs' else HERE/level
            dest.mkdir(exist_ok=True)
            metrics = extract(path, case, dest)
            for metric, value in metrics.items():
                metrics_rows.append((case, 'base' if level == 'runs' else level, metric, value))
            raw = HERE/'raw'/level/case
            raw.mkdir(parents=True, exist_ok=True)
            output_hashes = {}
            # Native numerical CSVs only: the executable also writes timing files.
            for src in sorted(path.glob(f'{case}_*.csv')):
                if '_summary_' not in src.name and '_profile_' not in src.name:
                    continue
                shutil.copyfile(src, raw/src.name)
                output_hashes[src.name] = digest(src)
            shutil.copyfile(deck, raw/deck.name)
            provenance['runs'].append(dict(case=case, level=level, command=[exe.name,deck.name],
                                          input_sha256=digest(deck), outputs_sha256=output_hashes,
                                          completed_successfully=True))
            print(case, level, metrics, flush=True)
    write_csv(HERE/'accuracy.csv', ['case','resolution','metric','error'], metrics_rows)
    (HERE/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    main()
