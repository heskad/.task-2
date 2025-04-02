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

## Контроль качества

### FastQC анализ
Перед запуском основного пайплайна необходимо выполнить контроль качества с помощью FastQC в два этапа:

1. **До тримминга и удаления адаптеров:**
```bash
fastqc -o qc/fastqc/raw/ data/raw/*.fastq
```

### Тримминг и удаление адаптеров
После первого FastQC анализа необходимо выполнить предварительную обработку чтений с помощью Trimmomatic:

```bash
trimmomatic PE -threads {threads} \
    data/raw/{sample}_R1.fastq data/raw/{sample}_R2.fastq \
    qc/trimmed/{sample}_R1.trimmed.fastq qc/trimmed/{sample}_R1.unpaired.fastq \
    qc/trimmed/{sample}_R2.trimmed.fastq qc/trimmed/{sample}_R2.unpaired.fastq \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
```

Параметры обработки:
- `ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10` - удаление адаптеров Illumina TruSeq3
- `SLIDINGWINDOW:4:20` - фильтрация по качеству с использованием скользящего окна
- `MINLEN:50` - удаление чтений короче 50 нуклеотидов

2. **После тримминга и удаления адаптеров:**
```bash
fastqc -o qc/fastqc/trimmed/ qc/trimmed/*.trimmed.fastq
```

Это позволит оценить качество данных на разных этапах обработки и убедиться в эффективности удаления адаптеров и фильтрации.

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
---

Snakemake -j <"кол-во задач"> 

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