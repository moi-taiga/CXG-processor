# CXG Processor Pipeline

A pipeline for pulling data from cellXgene and undertaking some basic preprocessing.

## User Story

As a single-cell Bioinformatitian, I want to automate my single-cell analysis, so that my workflow is faster.

## MVP

Automated pipeline to take CXG data and apply scanpy preprocessing.

## Stretch Goal

Continue the pipeline to GeneSwitches.

## Output Structure

Each run generates a directory containing:
- the data,
- a logfile,
- a config file,
- any graphical outputs (e.g. UMAPs, QC metrics),
- and the processed object saved as a `.h5ad` file.

## Structure

See `configs/` for configuration templates. Each run is stored in `analyses/`.

## Usage

1. Edit config files in `configs/`.
2. Run the pipeline via `src/main.py`.
3. Results are stored in a timestamped folder under `analyses/`.

## Requirements

See `environment.yml` or `requirements.txt`.
