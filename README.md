# UHBDB/toolkit
Unified Human Bacteriome Database

## Introduction
UHBDB/toolkit is a Nextflow pipeline for updating the Unified Human Bacteriome Database

## Databases

### Current database
| Release | URL | Size (GB) | Sequence count | Description | Samplesheet |
|---------|-----|-----------|----------------|-------------|-------------|
| 2026-03-11 | | 45 GB | 1,001,767 | HQ MAGs from mOTUs-db that were derived from human environments | |

### Database history
| Release | URL | Size (GB) | Sequence count | Description | Samplesheet |
|---------|-----|-----------|----------------|-------------|-------------|
| 2026-03-05 |  | 35 GB | 825,099 | Initial release of UHBDB created from Progenomes3, GTDB r214, and select JGI HQ isolate assemblies in human-associated genera. | |

## Pipeline summary
Below is a schematic overview of the UHVDB toolkit. For a more detailed explanation of each step, see the [`SUMMARY.md`](/.github/SUMMARY.md) file.

## Usage
> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

```bash
nextflow run UHVDB/toolkit -profile <singularity/institute>,test
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [`USAGE.md`](/.github/USAGE.md) and [`PARAMETERS.md`](/.github/PARAMETERS.md) files.

## Outputs

For details about the output files and reports, please refer to the [`OUTPUT.md`](/.github/OUTPUT.md) file.

## Future Goals

## Contributions and Support

## Credits

## Citations
