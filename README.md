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

## Arena

one-shot RefSeq/TOGA2-derived annotation of the dog genome (ROSCfam1.0):


| Tool                                                              | Threshold |       TP ↑ |       TN ↑ |   FP ↓ |    FN ↓ | Sensitivity ↑ | Specificity ↑ |       F1 ↑ | Balanced accuracy ↑ |  ROC AUC ↑ | Average precision ↑ |
| :---------------------------------------------------------------- | --------: | ---------: | ---------: | -----: | ------: | ------------: | ------------: | ---------: | ------------------: | ---------: | ------------------: |
| **[xORF](https://github.com/alejandrogzi/xorf)**                  |     0.500 | **24,947** |     16,335 |    118 | **360** |    **0.9858** |        0.9928 | **0.9905** |          **0.9893** | **0.9974** |          **0.9988** |
| [TranslationAI](https://github.com/rnasys/TranslationAI)          |     0.500 |     24,325 |     16,253 |    200 |     982 |        0.9612 |        0.9878 |     0.9763 |              0.9745 |     0.9878 |              0.9918 |
| [TD2 strict](https://github.com/Markusjsommer/TD2)                |     0.500 |     23,989 |     15,992 |    461 |   1,318 |        0.9479 |        0.9720 |     0.9642 |              0.9600 |     0.9716 |              0.9778 |
| [TRANSAID](https://github.com/wuzengding/TRANSAID)                |     0.635 |     23,401 | **16,354** | **99** |   1,906 |        0.9247 |    **0.9940** |     0.9589 |              0.9593 |     0.9861 |              0.9893 |
| [TD2](https://github.com/Markusjsommer/TD2)                       |     0.500 |     24,509 |     14,365 |  2,088 |     798 |        0.9685 |        0.8731 |     0.9444 |              0.9208 |     0.9793 |              0.9856 |
| [PSAURON](https://github.com/salzberg-lab/PSAURON)                |     0.500 |     24,844 |     13,688 |  2,765 |     463 |        0.9817 |        0.8319 |     0.9390 |              0.9068 |     0.9845 |              0.9908 |
| [TransDecoder](https://github.com/TransDecoder/TransDecoder)      |     0.500 |     17,986 |     13,394 |  3,059 |   7,321 |        0.7107 |        0.8141 |     0.7761 |              0.7624 |     0.7997 |              0.8161 |
| [GeneMark](http://topaz.gatech.edu/GeneMark/license_download.cgi) |     0.500 |     16,517 |     13,662 |  2,791 |   8,790 |        0.6527 |        0.8304 |     0.7404 |              0.7415 |     0.7723 |              0.8388 |


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
