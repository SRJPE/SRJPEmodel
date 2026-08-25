-- 0001_create_canonical_model_run.sql
--
-- Adds a table to track which model_run is "canonical" (the production fit)
-- for each model type + scope (site / run_year / site_selection, where
-- applicable). Append-only: never UPDATE or DELETE a row here — to change
-- the canonical run for a scope, INSERT a new row. The "current" designation
-- for a scope is whichever row has the latest set_at, mirroring the
-- newest-wins pattern already used for pin versions in model_run. This
-- gives a full audit trail for free, with no separate is_current flag or
-- transactional flip needed.
--
-- Note: no schema for model_name / model_run is currently tracked in git —
-- this repo's DB has been managed by hand against the live Postgres
-- instance. This file starts a db/migrations/ convention; it does not
-- attempt to also capture the pre-existing tables it references.
--
-- Run manually against the target database (staging/prod) — there is no
-- migration runner wired up yet.

CREATE TABLE canonical_model_run (
  id             SERIAL PRIMARY KEY,
  model_name_id  INTEGER NOT NULL REFERENCES model_name(id),
  site           TEXT,
  run_year       INTEGER,
  site_selection TEXT,
  model_run_id   INTEGER NOT NULL REFERENCES model_run(id),
  set_by         TEXT NOT NULL,
  set_at         TIMESTAMPTZ NOT NULL DEFAULT now(),
  note           TEXT
);

-- Speeds up "latest row per scope" queries (list_canonical_model_runs(),
-- get_canonical_history() in R/canonical.R).
CREATE INDEX canonical_model_run_scope_idx
  ON canonical_model_run (model_name_id, site, run_year, site_selection, set_at DESC);
