# Batch Processor for Step 2

This batch processor automates the processing of multiple simulation files from Step 1 through the Step 2 pipeline.

## Quick Start

```bash
# Process all simulations
python batch_processor.py

# Process specific rates
python batch_processor.py --rates 0.001,0.005,0.010

# Dry run (see what would be done)
python batch_processor.py --dry-run

# Resume from previous run
python batch_processor.py --resume
```

## Features

- **Automatic discovery** of simulation files
- **Rate-specific output directories** for organization
- **Resume capability** - pick up where you left off
- **Progress tracking** with ETA
- **Comprehensive logging** with timestamps
- **Error handling** - continues processing even if one simulation fails
- **State management** - tracks completed and failed simulations
- **Configurable** via YAML file

## Output Structure

```
step2/data/
├── rate_0.001000/
│   ├── snapshots/
│   │   └── year50_snapshot.json.gz
│   ├── lineages/
│   │   ├── mutant/
│   │   └── control/
│   └── plots/
│       └── year50_jsd_distribution.png
├── rate_0.002000/
│   └── ... (same structure)
└── .batch_state.json (processing state)
```

## Configuration

Edit `config/batch_config.yaml` to customize:
- Which year to extract
- Processing timeouts
- Logging settings
- Skip existing outputs
- And more...

## Logs

Check `logs/batch_run_YYYYMMDD_HHMMSS.log` for detailed processing information.

## Troubleshooting

- **Out of memory**: Reduce `parallel_jobs` in config (currently set to 1)
- **Script not found**: Ensure you're running from `step2/scripts/`
- **Permission denied**: Check write permissions on output directory

## Next Steps

After processing, run Step 3 to create mixed populations from these lineages.