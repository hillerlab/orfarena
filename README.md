<p align="center">
  <p align="center">
    <img width=100 align="center" src="./assets/figures/orfarena.png" >
  </p>

  <p align="center">
    <img width=200 align="center" src="./assets/figures/hillerlab.png" >
  </p>

  <span>
    <h1 align="center">
        orfarena
    </h1>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/orfarena" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/orfarena?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> benchmark ORF prediction methods </span>
        <br>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://en.wikipedia.org/wiki/Open_reading_frame">orf</a> .
        <a href="https://github.com/hillerlab/orfarena/blob/main/assets/pipeline/orfarena.mermaid">pipeline</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

---

## Usage

> [!NOTE]
> Requirements: Nextflow ≥ 25.04.6, Docker or Apptainer, Java.

```bash
git clone https://github.com/hillerlab/orfarena.git
cd orfarena
```

Edit `params.json` (set `genome`, `bed`), then:

```bash
# Docker
nextflow run main.nf -params-file params.json -profile docker

# Apptainer / Singularity
nextflow run main.nf -params-file params.json -profile apptainer
```

Smoke test:
```bash
nextflow run main.nf -profile test,apptainer
```

> [!NOTE]
> You can also specify these options directly in `params.json`.

---

## Output

```
results/
├── transdecoder/      *.bed
├── genemarkst/        *.bed
├── td2/               *.bed
├── td2_strict/        *.bed
├── psauron/           *.bed
├── gstama/            *.bed
└── pipeline_info/     timeline, trace, DAG
```

---

## Where to edit

| File | What |
|------|------|
| `params.json` | Paths, settings, checkpoints — per run |
| `nextflow.config` | Compute resources, profiles, container, SLURM — rarely |
