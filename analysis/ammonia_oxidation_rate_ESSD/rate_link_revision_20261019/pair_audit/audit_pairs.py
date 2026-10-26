"""Read-only source audit; writes only into this new audit directory.

Exact metadata matches are candidates, not proof of a shared sample/bottle.
Reconciliation with the cached original Tolar event table is separately flagged.
No regression, ecotype attribution, or future-rate prediction is performed.
"""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd

import os
if not os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the separate analysis workspace")
_workspace = Path(os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
OUT = _workspace / 'ammonia_oxidation_rate_ESSD/rate_link_revision_20261019/pair_audit'
OUT.mkdir(parents=True, exist_ok=True)
BASE = OUT.parents[1]
DB = BASE / 'Database of nitrification and nitrifiers in the global ocean.xlsx'
OLD = BASE / 'rate_link_mainline_20261003'
KEYS = ['study_key', 'date_text', 'obs_lat', 'obs_lon', 'sample_depth']


def save(frame, name):
    frame.to_csv(OUT/name, sep='\t', index=False, encoding='utf-8-sig')


def load(sheet):
    d = pd.read_excel(DB, sheet_name=sheet, engine='openpyxl')
    d.columns = d.columns.str.strip()
    d['excel_row'] = np.arange(len(d)) + 2
    raw = d['Data source'].fillna('').astype(str).str.strip()
    d['study'] = raw.where(raw.str.contains(r'et al|^Ward\b|^Ward\.|unpublished',
                                          regex=True, case=False)).ffill()
    d['study_key'] = d.study.str.strip().str.rstrip('. ').str.lower()
    d['date_text'] = d['Date'].fillna('').astype(str)
    mapping = {'Latitude':'obs_lat', 'Longitude':'obs_lon', 'Depth (m)':'sample_depth',
               'Rate (nmol N L-1 d-1)':'rate', 'total AOA (copies L-1)':'total_aoa'}
    for c in d.columns:
        for prefix, new in [('Temperature','temp'), ('Salinity','sal'),
                            ('Ammonium+ammonia','nh4'), ('O2','o2')]:
            if c.startswith(prefix):
                mapping[c] = new
    for old, new in mapping.items():
        if old in d:
            d[new] = pd.to_numeric(d[old], errors='coerce')
    d['obs_lon'] = (d.obs_lon+180)%360-180
    return d


def main():
    rates = load('volumetric ammonia oxidation')
    aoa = load('ammonia oxidizer')
    aoa = aoa[aoa.total_aoa.gt(0) & aoa.date_text.ne('')].dropna(subset=KEYS).copy()
    gridpath = BASE/'ecotype_function_20261002/paired_coastal_grids.tsv'
    grid = pd.read_csv(gridpath, sep='\t', usecols=['grid_id','Longitude','Latitude'])
    assert not grid.duplicated(['Longitude','Latitude']).any()
    sourcepath = OLD/'Tolar2020_BCODMO774848.csv'
    source = pd.read_csv(sourcepath, skiprows=[1])
    source['source_data_row'] = np.arange(len(source))+1
    sourcekeys = ['Cruise','Date','Station','depth']
    assert not source.duplicated(sourcekeys).any()
    for c in ['depth','Temperature','Salinity','NH4','Nit_rate_insitu','Nit_rate_dark',
              'amoA_qPCR','WCA_qPCR','WCB_qPCR','latitude','longitude']:
        source[c] = pd.to_numeric(source[c], errors='coerce')

    grouped = aoa.groupby(KEYS, as_index=False).agg(
        total_aoa_copies_L=('total_aoa','median'),
        aoa_rows=('excel_row','size'),
        aoa_excel_rows=('excel_row',lambda s:';'.join(map(str,s))),
        aoa_values=('total_aoa',lambda s:';'.join(map(str,s))),
        aoa_unique_values=('total_aoa','nunique'),
        qPCR_temp_C=('temp','median'))
    counts, studies, keyaudits, reconciled = [], [], [], []
    calibrationpath = BASE/'fixed_coefficient_scenario_20261016/source_included_100_records.tsv'
    calibration = pd.read_csv(calibrationpath, sep='\t')
    calibration_audit = None
    for maximum in [20,200]:
        r = rates[rates.rate.gt(0) & rates.sample_depth.between(0,maximum)
                  & rates.date_text.ne('')].dropna(subset=KEYS).copy()
        left_keys = r.groupby(KEYS).size().rename('rate_rows').reset_index()
        right_keys = aoa.groupby(KEYS).size().rename('aoa_rows').reset_index()
        k = left_keys.merge(right_keys, on=KEYS, validate='one_to_one')
        k['raw_join_rows'] = k.rate_rows*k.aoa_rows
        k['depth_limit_m'] = maximum
        keyaudits.append(k)
        paired = r.merge(grouped, on=KEYS, validate='many_to_one')
        paired['Longitude'] = np.floor((paired.obs_lon+180)/.25)*.25-180+.125
        paired['Latitude'] = np.floor((paired.obs_lat+90)/.25)*.25-90+.125
        paired = paired.merge(grid, on=['Longitude','Latitude'], how='left', validate='many_to_one')
        paired['date_is_exact_day'] = paired.date_text.str.match(r'^\d{4}-\d{2}-\d{2}(?: 00:00:00)?$')
        paired['same_bottle_verified'] = False
        paired['measured_LowT_HighT_available'] = False
        paired['original_source_event_verified'] = False
        paired['source_data_row'] = np.nan

        # Exact date and coordinate/depth agreement only locates the original event;
        # rates and amoA values must then both agree to flag source verification.
        target = paired[paired.study_key.eq('tolar et al., 2020. l&o')].copy()
        target['Date'] = target.date_text.str[:10]
        mapped = target.merge(source, on='Date', suffixes=('_tang','_source'))
        mapped = mapped[np.isclose(mapped.obs_lat,mapped.latitude)
                        & np.isclose(mapped.obs_lon,mapped.longitude)
                        & np.isclose(mapped.sample_depth,mapped.depth)].copy()
        assert not mapped.excel_row.duplicated().any()
        mapped['rate_verified'] = np.isclose(mapped.rate,mapped.Nit_rate_insitu,rtol=1e-6,atol=1e-8)
        mapped['amoA_verified'] = np.isclose(mapped.total_aoa_copies_L,mapped.amoA_qPCR,rtol=1e-6)
        verified = mapped[mapped.rate_verified & mapped.amoA_verified]
        lookup = dict(zip(verified.excel_row,verified.source_data_row_source))
        paired['source_data_row'] = paired.excel_row.map(lookup)
        paired['original_source_event_verified'] = paired.source_data_row.notna()
        mapped['depth_limit_m'] = maximum
        reconciled.append(mapped[['depth_limit_m','excel_row','aoa_excel_rows','grid_id',
                                 'source_data_row_source','Date','Cruise','Station','depth',
                                 'rate','Nit_rate_insitu','total_aoa_copies_L','amoA_qPCR',
                                 'rate_verified','amoA_verified']])
        keep = KEYS+['excel_row','rate','temp','sal','nh4','total_aoa_copies_L','aoa_rows',
                     'aoa_excel_rows','aoa_values','aoa_unique_values','qPCR_temp_C',
                     'Longitude','Latitude','grid_id','date_is_exact_day','same_bottle_verified',
                     'measured_LowT_HighT_available','original_source_event_verified','source_data_row']
        save(paired[keep], f'candidate_pairs_0_{maximum}m.tsv')
        if maximum == 20:
            coastal = paired[paired.grid_id.notna()]
            old = pd.read_csv(OLD/'coastal_rate_AOA_candidate_pairs_audited.tsv',sep='\t')
            check = coastal[['excel_row','rate','total_aoa_copies_L']].merge(
                old[['excel_row','rate','total_aoa_copies_L']], on='excel_row',
                how='outer', suffixes=('_new','_old'), indicator=True,
                validate='one_to_one')
            assert check['_merge'].eq('both').all()
            assert np.allclose(check.rate_new,check.rate_old)
            assert np.allclose(check.total_aoa_copies_L_new,check.total_aoa_copies_L_old)
            save(check, 'previous_43_candidate_reproduction_check.tsv')
            calibration_audit = calibration.merge(paired[['excel_row','total_aoa_copies_L',
                'aoa_excel_rows','original_source_event_verified']], on='excel_row',
                how='left', validate='one_to_one')
            save(calibration_audit,'fixed_coefficient_100_record_pair_status.tsv')
        for scope, p in [('all',paired),('strict_coastal_grid',paired[paired.grid_id.notna()])]:
            complete = p[p.date_is_exact_day & p.aoa_rows.eq(1) & p.temp.notna() & p.nh4.ge(0)]
            counts.append(dict(depth_limit_m=maximum,scope=scope,candidate_rate_rows=len(p),
                               event_keys=len(p[KEYS].drop_duplicates()),studies=p.study_key.nunique(),
                               coastal_grids=p.grid_id.nunique(),exact_day_rows=int(p.date_is_exact_day.sum()),
                               ambiguous_aoa_key_rows=int(p.aoa_rows.gt(1).sum()),
                               complete_temp_nh4_exact_unique=len(complete),
                               complete_studies=complete.study_key.nunique(),
                               complete_grids=complete.grid_id.nunique(),
                               source_event_verified=int(p.original_source_event_verified.sum()),
                               source_event_verified_studies=p.loc[p.original_source_event_verified,'study_key'].nunique(),
                               same_bottle_verified=0,measured_LowT_HighT_pairs=0))
            for study, q in p.groupby('study_key'):
                studies.append(dict(depth_limit_m=maximum,scope=scope,study=study,rows=len(q),
                                    event_keys=len(q[KEYS].drop_duplicates()),grids=q.grid_id.nunique(),
                                    exact_day_rows=int(q.date_is_exact_day.sum()),
                                    temperature_rows=int(q.temp.notna().sum()),
                                    ammonia_rows=int(q.nh4.ge(0).sum()),
                                    source_event_verified=int(q.original_source_event_verified.sum())))

    save(pd.DataFrame(counts),'pair_counts.tsv')
    save(pd.DataFrame(studies),'pair_studies.tsv')
    save(pd.concat(keyaudits,ignore_index=True),'join_key_multiplicity.tsv')
    save(pd.concat(reconciled,ignore_index=True),'tang_tolar_reconciliation.tsv')
    source_summary = []
    surface = source[source.depth.between(0,20)].copy()
    for name,c in [('in_situ_light','Nit_rate_insitu'),('dark','Nit_rate_dark')]:
        pairs = surface[surface[c].gt(0) & surface.amoA_qPCR.gt(0)]
        complete = pairs[pairs.NH4.ge(0)].dropna(subset=['Temperature','Salinity'])
        source_summary.append(dict(assay=name,positive_rate_amoA_pairs=len(pairs),
                                   complete_environment_rows=len(complete),
                                   cruises=complete.Cruise.nunique(),stations=complete.Station.nunique(),
                                   source_event_key_duplicates=int(complete.duplicated(sourcekeys).sum()),
                                   measured_LowT_HighT_pairs=0,same_bottle_verified=0))
        save(complete, f'original_tolar_surface_{name}_complete_rows.tsv')
    q = source.dropna(subset=['amoA_qPCR','WCA_qPCR','WCB_qPCR'])
    q = q[q.amoA_qPCR.gt(0)]
    excluded = q[q.WCA_qPCR+q.WCB_qPCR>q.amoA_qPCR]
    save(excluded[sourcekeys+['amoA_qPCR','WCA_qPCR','WCB_qPCR']], 'qPCR_nonadditivity.tsv')
    audit = pd.concat(keyaudits, ignore_index=True)
    metadata_rows = []
    for sheet in ['rate metadata','qPCR metadata']:
        metadata = pd.read_excel(DB,sheet_name=sheet,engine='openpyxl')
        for idx,row in metadata.iterrows():
            study = str(row.iloc[0])
            if 'Tolar' in study and '2020' in study:
                for field,value in row.items():
                    if pd.notna(value):
                        metadata_rows.append(dict(sheet=sheet,excel_row=idx+2,
                                                  study=study,field=field,value=value))
    save(pd.DataFrame(metadata_rows),'tolar_workbook_assay_metadata.tsv')
    multiplicity = []
    for limit,k in audit.groupby('depth_limit_m'):
        multiplicity.append(dict(depth_limit_m=int(limit),matched_keys=len(k),
                                 repeated_rate_keys=int(k.rate_rows.gt(1).sum()),
                                 repeated_aoa_keys=int(k.aoa_rows.gt(1).sum()),
                                 many_to_many_keys=int((k.rate_rows.gt(1)&k.aoa_rows.gt(1)).sum()),
                                 raw_join_rows=int(k.raw_join_rows.sum()),
                                 median_aoa_join_rows=int(k.rate_rows.sum())))
    summary = dict(date='2026-10-19',primary_workbook=str(DB),pair_counts=counts,
                   multiplicity=multiplicity,original_tolar_rows=len(source),
                   original_tolar_surface_rows=len(surface),original_tolar_source_summary=source_summary,
                   qPCR_all_three_assays=len(q),qPCR_nonadditive=len(excluded),
                   fixed_coefficient_calibration=dict(records=len(calibration),
                       min_sample_depth_m=float(calibration.sample_depth.min()),
                       max_sample_depth_m=float(calibration.sample_depth.max()),
                       measured_total_AOA_candidate_records=int(calibration_audit.total_aoa_copies_L.notna().sum()),
                       source_event_verified=int(calibration_audit.original_source_event_verified.eq(True).sum())),
                   decision='No future-rate model: no observed LowT/HighT rate pairs and no compatible future total-amoA inputs.',
                   evidence_limits=['Metadata keys define candidate events, not bottle identity.',
                                    'Original Tolar table verifies co-recorded event measurements, not same aliquot.',
                                    'Total archaeal amoA copies/L is not ecotype abundance or measured cells/L.',
                                    'Surface total-AOA pilot is regional; it cannot calibrate future relative LowT/HighT predictions.',
                                    'Positive-rate and positive-abundance eligibility excludes zero/BDL observations.'])
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    files=[DB,gridpath,sourcepath,calibrationpath,Path(__file__).resolve()]
    save(pd.DataFrame([dict(path=str(p),bytes=p.stat().st_size,
                           sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]), 'input_manifest.tsv')
    print(json.dumps(summary,indent=2,ensure_ascii=False))


if __name__=='__main__':
    main()
