#  Метагеномный анализ эукариот (на примере миноги)

**Snakemake-пайплайн для сборки, бининга и аннотации метагеномных данных**

---

## Обзор
Автоматизированный пайплайн для анализа метагеномных данных эукариот (например, миноги *Petromyzon marinus*) с поддержкой:
- Illumina и Nanopore данных
- Контроля качества
- *De novo* сборки
- Таксономической классификации
- Аннотации генов

---

## Установка

### Предварительные требования
- [Miniconda/Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
- Рекомендуется Mamba для ускорения:

```bash
conda install -n base -c conda-forge mamba
---

git clone https://github.com/heskad/.task-2
cd .task-2
cd envs
chmod +x install_envs.sh
./install_envs.sh
```


## Этапы анализа

| Этап анализа          | Инструменты                          | Выходные файлы                     |
|-----------------------|--------------------------------------|------------------------------------|
| **Контроль качества** | FastQC, Trimmomatic, Cutadapt        | `results/qc/trimmed/*.fastq`       |
| **Сборка генома**     | MEGAHIT (Illumina),             | `results/assembly/final_assembly.fasta` |
| **Таксономия**        | Kraken2, Kaiju, BUSCO                | `results/taxonomy/kraken_out.txt`  |
| **Биннинг**           | MetaBAT2, CONCOCT                    | `results/binning/metabat2/bins/`   |
| **Аннотация**         | Prokka, BLAST                        | `results/annotation/prokka.gff`    |
| **Визуализация**      | QUAST, MultiQC                       | `results/quast/report.html`        |