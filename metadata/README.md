# Sample metadata

`metadata/sample_metadata.tsv` is generated from the public GEO family SOFT file rather than hand-entered.

The parser preserves GEO fields including:

- `sample_id`
- `title`
- `source_name`
- `nas_score`
- `fibrosis_stage`
- `group_in_paper`
- `disease`
- `stage`
- `description`

It also derives `analysis_group` for a transparent fibrosis-focused comparison:

- `Control` for GEO samples labeled as controls;
- `NAFLD_F0_F1` for NAFLD samples with fibrosis stage 0–1;
- `NAFLD_F2` for NAFLD samples with fibrosis stage 2;
- `NAFLD_F3_F4` for NAFLD samples with fibrosis stage 3–4.

The derived grouping does **not** replace or alter the original GEO metadata. It exists only as an analysis convenience and is created deterministically from the public `disease` and `fibrosis_stage` fields.
